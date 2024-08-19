using Interpolations
using CategoricalArrays, DataFrames, XLSX
using Plots
using ProgressMeter
using Statistics
using DSP
using Impute
using CubicSplines
using Feather


function calculate_sr(df::DataFrame)
    # calculate the sampling rate by getting the ms difference between the first and last sample of each event/subject/eye and averaging the difference
    df_grouped = groupby(df, [:subject, :eye, :event_id])
    sr = []
    for (key, subdf) in pairs(df_grouped)
        sr_subj = 1000 / mean(diff(subdf.ms))
        push!(sr, sr_subj)
    end
    return mean(sr)
end


function upsample(df::DataFrame, target_sr::Int, group_cols::Vector{Symbol})
    df_grouped = groupby(df, group_cols)
    df_upsampled = DataFrame()
    @showprogress "Upsampling to $target_sr Hz" for (key, subdf) in pairs(df_grouped)
        t = subdf.ms
        v = subdf.diameter_px
        aud = subdf.audio
        blink = subdf.blink
        conf = subdf.confidence
        artifact = subdf.artifact
        itp = linear_interpolation(t, v)
        itp_aud = constant_interpolation(t, aud)
        itp_blink = constant_interpolation(t, blink)
        itp_conf = linear_interpolation(t, conf)
        itp_artifact = constant_interpolation(t, artifact)
        t_start = t[begin]
        t_end = t[end]
        t_new = t_start:1000/target_sr:t_end
        v_new = itp(t_new)
        v_new_len = length(v_new)
        aud_new = itp_aud(t_new)
        blink_new = itp_blink(t_new)
        conf_new = itp_conf(t_new)
        artifact_new = itp_artifact(t_new)
        df_upsampled = vcat(df_upsampled, DataFrame(
            subject = repeat([key.subject], v_new_len),
            eye = repeat([key.eye], v_new_len),
            event_id = repeat([key.event_id], v_new_len),
            ms = t_new,
            diameter_px = v_new,
            audio = aud_new,
            blink = blink_new,
            confidence = conf_new,
            artifact = artifact_new
        ))

    end
    return df_upsampled
end

function deblink_smooth_interpolate(df::DataFrame; sr::Int = 1000, extend::Int = 16, window_length::Int = 32)
    df_grouped = groupby(df, [:subject, :eye, :event_id])
    df_deblinked = DataFrame()
    if mod(window_length, 2) == 0
        window_length += 1
    end
    @showprogress "Deblinking, etc." for (key, subdf) in pairs(df_grouped)
        diam = Vector{Union{Missing, Float64}}(subdf.diameter_px)
        indices = (subdf.blink .== 1) .| (subdf.artifact .== 1)
        if any(indices)

            if extend > 0
                for i in 1:extend
                    
                    indices = indices .| 
                        circshift(
                            vcat(
                                indices[begin:end-i],
                                fill(false, i)
                            ),
                            i) .|
                        circshift(
                            vcat(
                                fill(false, i),
                                indices[1+i:end]
                            ),
                            -i
                        )
                end
                indices[begin] = false
                indices[end] = false
            end
            diam[indices] .= missing
            diam = Impute.interp(diam) |> Impute.locf() |> Impute.nocb()
        end

        # smoothing
        w = hanning(window_length)
        y = Float64.(diam)
        y_new = Float64.([2*y[1].-reverse(y[1:window_length]); y[:]; 2*y[end].-reverse(y[end-window_length:end])])
        y_smooth = conv(y_new, w/sum(w))
        ind = floor(Int, 1.5*window_length)
        y_smooth = y_smooth[1+ind:end-ind-1]

        # cubic spline Interpolation
        if any(indices)
            x = Float64.(subdf.ms[.!indices])
            y = Float64.(y_smooth[.!indices])
            try 
                cs = CubicSpline(x, y)
                y_smooth[indices] .= cs(Float64.(subdf.ms[indices]))
            catch e
                println("Error in CubicSpline for $key, rejecting")
                continue
            end
            
        end

        first_audio_idx = findfirst(subdf.audio .== 1)
        df_deblinked = vcat(df_deblinked, DataFrame(
            subject = repeat([key.subject], length(y_smooth)),
            eye = repeat([key.eye], length(y_smooth)),
            event_id = repeat([key.event_id], length(y_smooth)),
            ms = subdf.ms,
            relative_ms = collect(
                -(first_audio_idx - 1) * 1000/sr:1000/sr:(nrow(subdf) - first_audio_idx) * 1000/sr
            ),
            diameter_px = y_smooth,
            audio = subdf.audio,
            blink = subdf.blink,
            confidence = subdf.confidence,
            artifact = subdf.artifact,
            extend_mask = indices,
        ))
    end
    return df_deblinked
end
