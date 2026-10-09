# Production passage inference — reference-only v3

This folder implements the final requested Bridge West workflow:

**A031/B053 from each of 12 passages in 2022 and 24 passages in 2025 -> infer passage loads and condition the frozen shaker-informed Stabil discrepancy -> predict A058/A070/B031/B081.**

Target measurements are used only after prediction for assessment.

## Inputs

Place these original archives below one data root; subfolders are allowed:

- `New folder (2).zip` — metadata
- `ltw21_1.zip`, `ltw21_2.zip`
- `ltw24_1.zip`, `ltw24_2.zip`, `ltw24_3.zip`, `ltw24_4.zip`

The ST-W ZIPs are not required at prediction time because the frozen 84-state shaker ensemble is included in the package.

## Run all 36 passages

Edit the two paths in `matlab/run_all_bridgewest_stabil_v3.m`, then run it, or call the stages explicitly:

```matlab
run_bridgewest_stabil_shaker_informed(dataRoot,outputRoot,'Campaign','all');
score_bridgewest_stabil_results(dataRoot,outputRoot);
plot_bridgewest_stabil_results(outputRoot,'Event','2025_run40');
```

MATLAB needs Signal Processing Toolbox for the 10 Hz zero-phase strain preprocessing.

## Inference/assessment firewall

Prediction stage:

- reads bridge strain A031 and B053 only;
- obtains speed/entry time from the independent rail timing channels;
- infers eight bounded bogie resultants;
- conditions the 84-state shaker-informed structural discrepancy;
- predicts A058, A070, B031 and B081;
- saves the prediction before any target measurement is opened.

Assessment stage:

- opens A058/A070/B031/B081 only after the prediction file exists;
- reuses the already-frozen speed and entry time;
- cannot modify loads, model parameters, discrepancy states or prediction intervals;
- writes both all-target and primary-quality score tables.

The production protocol does **not** fit a second target-derived/local discrepancy model.

## Main output folders

- `processed_references/` — preprocessing records containing only A031/B053;
- `permitted_inputs/` — exact two-channel inference input for each passage;
- `predictions/` — all four target predictions plus load posterior and displacement predictions;
- `validation_targets/` — created only by the post-prediction scoring stage;
- `tables/target_assessment_all.csv` — all four target rows for every predicted passage;
- `tables/target_scores_primary.csv` — only quality-`ok` target observations;
- `tables/target_summary_primary.csv` — campaign/target aggregate of the primary rows;
- `tables/target_summary_uninspected_only.csv` — sensitivity summary excluding previously inspected run 40;
- `tables/reference_scores.csv` — A031/B053 posterior-fit diagnostics.

Run 40 remains included in the requested 24-passage 2025 set and is explicitly flagged as previously inspected.
