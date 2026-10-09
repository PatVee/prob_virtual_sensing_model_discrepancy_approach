# MATLAB production script structure — v3

`run_bridgewest_stabil_shaker_informed.m` is the inference runner for the default 36 passages. For each passage it reads independent rail timing information, then preprocesses **only A031 and B053** as structural measurements. The target channels A058/A070/B031/B081 are never read by this stage.

`bw_predict_continuous_shaker.m` performs the fixed-point reference calculation. It projects the 84-state unit-bogie shaker discrepancy basis to the current eight-bogie mean load, computes discrepancy mean/cross-covariance, obtains a bounded continuous load posterior from A031/B053, and propagates the final load/discrepancy uncertainty to every model output.

`bw_continuous_load_posterior.m`, `bw_bounded_gaussian_moments.m` and `bw_truncnorm_scalar.m` implement the bounded eight-bogie posterior. The hard load bounds are 180–360 kN. No target measurement or 2025 weighing value enters this calculation.

`raw_functions/bw_prepare_raw_event.m` is deliberately strict: it requires an explicit list of permitted inference channels. The production runner passes only the two configured references. It saves `processed_references/<event>.mat` and `permitted_inputs/<event>.mat`, both of which contain A031/B053 only.

`score_bridgewest_stabil_results.m` is a separate post-prediction program. It calls `raw_functions/bw_prepare_validation_targets.m`, which reopens the raw passage only after the prediction exists and extracts A058/A070/B031/B081 using the already-frozen speed and entry time. All four target predictions receive assessment rows when measurement data are numerically available. Only quality-`ok` rows enter the primary aggregate; this preserves 2022 A070 predictions while respecting their reduced-quality flag.

`audit_bridgewest_reference_only_protocol.m` verifies the configured 12+24 event set, exactly two saved inference channels, exact reference names, absence of target names in permitted inputs, and saved prediction flags confirming no target use.

`plot_bridgewest_stabil_results.m` reads saved score tables and post-prediction target files only. It never calls inference or raw target extraction.

The detailed mechanics remain under package-level `matlab/` and `shaker_informed/matlab/`. Surface strain is obtained from the Stabil shell observation operator; the old equivalent-beam fitted strain-recovery depth is not used.
