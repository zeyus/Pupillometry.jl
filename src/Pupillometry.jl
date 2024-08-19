include("util/Data.jl")
include("util/Processing.jl")

all_subj_fixed = load_data("G:/data/pupillometry", use_cached=false, cache_file="G:/github_projects/Pupillometry/data/audio_segments_upscaled_smoothed.feather")

# plot the data
function mean_baseline_subtract(x; sr::Int = 1000)
    # find the first audio event
    fa = findfirst(x.audio .== 1)
    baseline = [
        x.diameter_px[idx] for (idx,ts) in enumerate(x.ms) if idx < fa
    ]
    if isempty(baseline)
        return x.diameter_px
    end
    return x.diameter_px .- mean(baseline)
end
# subtract the mean for time < 0 for each subject and eye from the diameter_px
audio_segments_df_grouped_subj = groupby(all_subj_fixed, [:subject, :eye, :event_id])
audio_segments_df_grouped_subj = combine(
    audio_segments_df_grouped_subj,
    AsTable([
        :diameter_px,
        :ms,
        :audio,
    ]) => (
        (x) -> mean_baseline_subtract(x)
    ) => :relative_diameter_px_scaled,
    :color,
    :eye,
    :ms,
    :ms => eachindex => :relative_ms, # only works at SR = 1000
    :diameter_px,
    # :diameter_px_orig,
    :event_id,
    :subject,
    :confidence,
    :artifact,
    :blink,
)

min_ms = 500
max_ms = 5000

audio_segments_df_grouped_subj = audio_segments_df_grouped_subj[
    audio_segments_df_grouped_subj.relative_ms .>= min_ms,
    :
]

audio_segments_df_grouped_subj = audio_segments_df_grouped_subj[
    audio_segments_df_grouped_subj.relative_ms .<= max_ms,
    :
]

# plot all events for one subject
example_subject = "2F"
single_subject_events = audio_segments_df_grouped_subj[audio_segments_df_grouped_subj.subject .== example_subject, :]
single_subject_events_left = single_subject_events[single_subject_events.eye .== "left", :]
single_subject_events_right = single_subject_events[single_subject_events.eye .== "right", :]
p1 = plot(
    single_subject_events_left.relative_ms,
    single_subject_events_left.relative_diameter_px_scaled,
    group = single_subject_events_left.event_id,
    alpha=0.5,
    legend=:none,
    xlabel="Time (ms)",
    ylabel="Pupil diameter (range scaled)",
    title="Single Subject Pupil diameter (L)"
)
# add artifact and blink markers
single_subject_blinks = single_subject_events_left[single_subject_events_left.blink .== 1, :]
single_subject_artifacts = single_subject_events_left[single_subject_events_left.artifact .== 1, :]
scatter!(single_subject_blinks.relative_ms, single_subject_blinks.relative_diameter_px_scaled, label="Blink", color=:red, markerstrokewidth=0.2, markerstrokealpha=0.1)
scatter!(single_subject_artifacts.relative_ms, single_subject_artifacts.relative_diameter_px_scaled, label="Artifact", color=:green, markerstrokewidth=0.2, markerstrokealpha=0.1)

p2 = plot(
    single_subject_events_right.relative_ms,
    single_subject_events_right.relative_diameter_px_scaled,
    group = single_subject_events_right.event_id,
    alpha=0.5,
    legend=:none,
    xlabel="Time (ms)",
    ylabel="Pupil diameter (range scaled)",
    title="Single Subject Pupil diameter (R)"
)
# add artifact and blink markers
single_subject_blinks = single_subject_events_right[single_subject_events_right.blink .== 1, :]
single_subject_artifacts = single_subject_events_right[single_subject_events_right.artifact .== 1, :]
scatter!(single_subject_blinks.relative_ms, single_subject_blinks.relative_diameter_px_scaled, label="Blink", color=:red, markerstrokewidth=0.2, markerstrokealpha=0.1)
scatter!(single_subject_artifacts.relative_ms, single_subject_artifacts.relative_diameter_px_scaled, label="Artifact", color=:green, markerstrokewidth=0.2, markerstrokealpha=0.1)



plot(p1, p2, layout=(2,1), size=(2200, 2200), thickness_scaling=2)

# save
savefig("single_subject_pupil_diameter.png")

audio_segments_df_grouped = groupby(audio_segments_df_grouped_subj, [:subject, :eye, :relative_ms])
audio_segments_df_grouped = combine(
    audio_segments_df_grouped,
    :relative_diameter_px_scaled => mean => :diameter_px_mean,
    :relative_diameter_px_scaled => std => :diameter_px_std,
    :color => first => :color,
    :eye => first => :eye,
    :relative_ms => first => :relative_ms,
    :subject => first => :subject
)

# filter NaNs
audio_segments_df_grouped = audio_segments_df_grouped[.!isnan.(audio_segments_df_grouped.diameter_px_std), :]

# crop to max(min(relative_ms)) and min(max(relative_ms)) where there is data for both eyes

min_l = minimum(audio_segments_df_grouped.relative_ms[audio_segments_df_grouped.eye .== "left"])
max_l = maximum(audio_segments_df_grouped.relative_ms[audio_segments_df_grouped.eye .== "left"])
min_r = minimum(audio_segments_df_grouped.relative_ms[audio_segments_df_grouped.eye .== "right"])
max_r = maximum(audio_segments_df_grouped.relative_ms[audio_segments_df_grouped.eye .== "right"])

min_ms = max(min_l, min_r)
max_ms = min(max_l, max_r)
max_ms = 5000

audio_segments_df_grouped = audio_segments_df_grouped[
    audio_segments_df_grouped.relative_ms .>= min_ms,
    :
]

audio_segments_df_grouped = audio_segments_df_grouped[
    audio_segments_df_grouped.relative_ms .<= max_ms,
    :
]

audio_segments_df_grouped_left = audio_segments_df_grouped[audio_segments_df_grouped.eye .== "left", :]
audio_segments_df_grouped_right = audio_segments_df_grouped[audio_segments_df_grouped.eye .== "right", :]

p1 = plot(
    audio_segments_df_grouped_left.relative_ms,
    audio_segments_df_grouped_left.diameter_px_mean,
    ribbon=audio_segments_df_grouped_left.diameter_px_std,
    group = audio_segments_df_grouped_left.subject,
    alpha=0.5,
    legend = :bottomleft,
    legendcolumns=5,
    fillalpha=0.1,
    xlabel="Time (ms)",
    ylabel="Pupil diameter (range scaled)",
    title="Pupil diameter over time (Left eye)"
)
p1 = hline(p1, [0.0], color=:black, alpha=0.5, linestyle=:dash, label="Baseline")
# add line for audio onset
p1 = vline(p1, [1000.0], color=:blue, alpha=0.5, linestyle=:dash, label="Audio onset")
p2 = plot(
    audio_segments_df_grouped_right.relative_ms,
    audio_segments_df_grouped_right.diameter_px_mean,
    ribbon=audio_segments_df_grouped_right.diameter_px_std,
    group = audio_segments_df_grouped_right.subject,
    alpha=0.5,
    legend=:none,
    fillalpha=0.10,
    xlabel="Time (ms)",
    ylabel="Pupil diameter (range scaled)",
    title="Pupil diameter over time (Right eye)"
)
# add baseline line
p2 = hline(p2, [0.0], color=:black, alpha=0.5, linestyle=:dash, label="Baseline")
p2 = vline(p2, [1000.0], color=:blue, alpha=0.5, linestyle=:dash, label="Audio onset")
plot(p1, p2, layout=(2,1), size=(2200, 2200), thickness_scaling=2)

# save plot 
savefig("pupil_diameter.png")

subj_10f = audio_segments_df_grouped[audio_segments_df_grouped.subject .== "10F", :]

sort!(subj_10f, :relative_ms)
# find how many values of relative_ms have more than 2 values
# count the number of unique values of relative_ms
n_unique = length(unique(subj_10f.relative_ms))
n_total = nrow(subj_10f)
n_total/n_unique


# plot 
plot(
    subj_10f.relative_ms,
    subj_10f.diameter_px_mean,
    ribbon=subj_10f.diameter_px_std,
    group = subj_10f.color,
    alpha=0.5,
    fillalpha=0.1,
    legend=:topleft,
    xlabel="Time (ms)",
    ylabel="Pupil diameter (px)",
    title="Pupil diameter over time"
)



audio_segments_df_grouped_subj_left = audio_segments_df_grouped_subj[audio_segments_df_grouped_subj.eye .== "left", :]

# plot the data for each subject left eye
plot(
    audio_segments_df_grouped_subj_left.relative_ms,
    audio_segments_df_grouped_subj_left.relative_diameter_px,
    group = audio_segments_df_grouped_subj_left.subject,
    alpha=0.5,
    legend=:topleft,
    xlabel="Time (ms)",
    ylabel="Pupil diameter (px)",
    title="Pupil diameter over time"
)
