using Interpolations
using CategoricalArrays, DataFrames, XLSX
using Plots
using ProgressMeter
using Statistics
using DSP
using Impute
using CubicSplines
using Feather
using Logging

function read_pupil_data(file::String, sheet::String = "audio")
    df = DataFrame(XLSX.readtable(file, sheet))
    return df
end

"""
    collate_subject_data(subject::String, basedir::String)

    Subjects are in a folder with the subject ID + their gender,
    e.g. 
      - 1M
      - 1F
      - 2M
      - 2F
      - ...
    Each folder contains two Excel files, one for the left eye (_dx) and one for the right eye (_sx).
    The data is in an Excel file with the following sheets:
    "audio" -> the audio condition
    "baseline" -> the baseline condition
"""
function collate_subject_data(subject::String, basedir::String)
    left_eye_audio = read_pupil_data(joinpath(basedir, subject, "$(subject)_dx.xlsx"))
    right_eye_audio = read_pupil_data(joinpath(basedir, subject, "$(subject)_sx.xlsx"))
    left_eye_baseline = read_pupil_data(joinpath(basedir, subject, "$(subject)_dx.xlsx"), "baseline")
    right_eye_baseline = read_pupil_data(joinpath(basedir, subject, "$(subject)_sx.xlsx"), "baseline")
    # rename columns
    col_rename_map = Dict(
        :"diameter (pixel)" => :"diameter_px",
        :"blink (1 yes/0 no)" => :blink,
        :"blink (1 sì/0 no)" => :blink,
        :"audio (1 yes/0 no)" => :audio,
        :"audio (1 sì/0 no)" => :audio,
        :"artifact (1 yes/0 no)" => :artifact,
        :"artefatto (1 sì/0 no)" => :artifact,
    )
    @showprogress for ds in [left_eye_audio, right_eye_audio, left_eye_baseline, right_eye_baseline]
        for (old_col, new_col) in col_rename_map
            if old_col in names(ds)
                rename!(ds, old_col => new_col)
            end
        end
    end

    return left_eye_audio, right_eye_audio, left_eye_baseline, right_eye_baseline
end

"""
    collate_all_subject_data(basedir::String)

    Collate all the data for all the subjects in the `basedir` directory.
    Return a single dataframe with the data and a column for subject, eye, and condition.
"""
function collate_all_subject_data(basedir::String)
    subjects = readdir(basedir)
    all_data = DataFrame()
    @showprogress for subject in subjects
        if !isdir(joinpath(basedir, subject))
            continue
        end
        # subjects 1M/F to 5M/F are LL lighting, 6M/F to 10M/F are HL lighting
        lighting = "HL"
        if parse(Int, rstrip(subject, ['M', 'F'])) <= 5
            lighting = "LL"
        end
        left_eye_audio, right_eye_audio, left_eye_baseline, right_eye_baseline = collate_subject_data(subject, basedir)
        left_eye_audio.subject .= subject
        left_eye_audio.eye .= "left"
        left_eye_audio.condition .= "audio"
        left_eye_audio.lighting .= lighting

        right_eye_audio.subject .= subject
        right_eye_audio.eye .= "right"
        right_eye_audio.condition .= "audio"
        right_eye_audio.lighting .= lighting

        left_eye_baseline.subject .= subject
        left_eye_baseline.eye .= "left"
        left_eye_baseline.condition .= "baseline"
        left_eye_baseline.lighting .= lighting

        right_eye_baseline.subject .= subject
        right_eye_baseline.eye .= "right"
        right_eye_baseline.condition .= "baseline"
        right_eye_baseline.lighting .= lighting

        all_data = vcat(
            all_data,
            left_eye_audio,
            right_eye_audio,
            left_eye_baseline,
            right_eye_baseline
        )
    end
    return DataFrame(all_data)
end


# audio column = 0 means no audio, 1 means audio played, there are multiple audio events per subject
# segment the data into ms_before_audio + duration of audio (while audio = 1) + ms_after_audio
# dataframe is a grouped dataframe by subject
function segment_data(df::GroupedDataFrame, ms_before_audio::Int, ms_after_audio::Int, sr::Int = 62)
    audio_segments::Dict{String, Vector{Any}} = Dict(
        "subject" => String[],
        "diameter_px" => AbstractFloat[],
        "diameter_px_orig" => AbstractFloat[],
        "eye" => String[],
        "audio" => Int[],
        "event_id" => Int[],
        "ms" => Int[],
        "relative_ms" => Int[],
        "confidence" => Float64[],
        "artifact" => Int[],
        "blink" => Int[],
    )
    n_subj = length(keys(df))
    p1 = Progress(n_subj; showspeed = true, desc = "Subject", color=:red, offset=2)
    start_frame_offset = floor(Int, ms_before_audio/1000 * sr)
    end_frame_offset = ceil(Int, ms_after_audio/1000 * sr)
    for (key, subdf) in pairs(df)
        # first calculate subject diameter dynamic range
        e_min = minimum(subdf.diameter_px)
        e_max = maximum(subdf.diameter_px)
        
        # scale the diameter to the range of the subject
        subdf.diameter_px_scaled = map(
            subdf.diameter_px,
        ) do diameter
            return (diameter - e_min) / (e_max - e_min)
        end

        audio_events = findall(subdf.audio .== 1)
        if isempty(audio_events)
            next!(p1)
            continue
        end
        # remove seuqential audio events
        firsts = diff(audio_events) .> 1
        firsts = vcat(BitVector([1]), firsts)
        audio_events = audio_events[firsts]
        n_events = length(audio_events)
        p2 = Progress(n_events; showspeed = true, desc = "Event", color=:blue, offset=1)
        event_id = 1
        
        for event in audio_events
            start_idx = max(1, event - start_frame_offset)
            end_idx = min(nrow(subdf), event + end_frame_offset)
            relative_ms = collect(-(event-start_idx)*1000/sr:1000/sr:(end_idx-event)*1000/sr)
            
            append!(audio_segments["subject"], repeat([key.subject], end_idx - start_idx + 1))
            append!(audio_segments["eye"], repeat([key.eye], end_idx - start_idx + 1))
            append!(audio_segments["diameter_px_orig"], subdf[start_idx:end_idx, :diameter_px])
            append!(audio_segments["diameter_px"], subdf[start_idx:end_idx, :diameter_px_scaled])
            append!(audio_segments["audio"], subdf[start_idx:end_idx, :audio])      
            append!(audio_segments["event_id"], repeat([event_id], end_idx - start_idx + 1))
            append!(audio_segments["ms"], subdf[start_idx:end_idx, :milliseconds])
            append!(audio_segments["relative_ms"], relative_ms)
            append!(audio_segments["confidence"], subdf[start_idx:end_idx, :confidence])
            append!(audio_segments["artifact"], subdf[start_idx:end_idx, :artifact])
            append!(audio_segments["blink"], subdf[start_idx:end_idx, :blink])

            event_id += 1
            next!(p2)
        end
        next!(p1)
    end
    finish!(p1)
    return DataFrame(audio_segments)
end

function update_columns(df::DataFrame)
    # make subject and eye categorical
    df.subject = CategoricalArray(df.subject)
    df.eye = CategoricalArray(df.eye)
    df.color = map(eye -> eye == "left" ? :blue : :red, df.eye)
    # make blink, artifact and audio boolean
    df.blink = map(blink -> blink == 1 ? true : false, df.blink)
    df.artifact = map(artifact -> artifact == 1 ? true : false, df.artifact)
    df.audio = map(audio -> audio == 1 ? true : false, df.audio)
    return df
end

function load_data(path::String; use_cached::Bool = true, lighting = "LL",  cache_file::String = "./data/audio_segments_upscaled_smoothed.feather")
    @info "Loading data from $path"
    @info "Looking for cache file $cache_file"
    if use_cached && isfile(cache_file)
        @info "Using cached data"
        all_subj_fixed = Feather.read(cache_file)
        return update_columns(all_subj_fixed)
    end
    @info "Cache file not found, processing data"
    @info "Collating all subject data"
    pupil_data = collate_all_subject_data(path)
    audio_data_LL = pupil_data[
        (pupil_data.condition .== "audio") .& (pupil_data.lighting .== lighting), :]

    audio_data = groupby(audio_data_LL, [:subject, :eye])

    ms_before_audio = 1000 
    ms_after_audio = 5000
    audio_segments_df = segment_data(audio_data, ms_before_audio, ms_after_audio)
    audio_segments_df_upsampled = upsample(audio_segments_df, 1000, [:subject, :eye, :event_id])
    all_subj_fixed = deblink_smooth_interpolate(
        audio_segments_df_upsampled;
        extend=100,
        window_length=65,
        sr=1000,
    )
    # now is a good time to save
    Feather.write(cache_file, all_subj_fixed)
    return update_columns(all_subj_fixed)
end
