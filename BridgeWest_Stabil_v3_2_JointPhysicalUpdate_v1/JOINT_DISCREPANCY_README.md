# Bridge West v3.2 — joint load/strain/displacement discrepancy branch

This branch addresses the amplitude-identifiability issue in the original shaker-only v3.2 model.

## Information split

- **2024 shaker tests:** frequency, damping, and force-normalised modal-residue evidence. The residue correction changes global force-to-response flexibility and therefore affects both strain and displacement.
- **2022 traffic campaign:** independent static bogie weighing plus **A031/B053 only** calibrate one scalar strain-recovery factor. No 2022 target strain is used.
- **2025 traffic campaign:** A031/B053 infer passage-specific bogie loads. 2025 static weighing, A058/A070/B031/B081, and PLS displacement are opened only after all predictions are frozen.

## Discrepancy layers

The joint covariance is the sum of four independent components:

1. 84 shaker frequency/damping states;
2. three distributed two-rail load-transfer footprints (`hx = 0.6, 1.2, 2.4 m`, `hy = 0.30, 0.45, 0.60 m`), equal epistemic weights;
3. common response/compliance gains `0.8, 1.0, 1.2` with weights `0.25, 0.50, 0.25`;
4. residual strain-recovery gains `0.92, 1.00, 1.08` with the same weights.

The added compliance and recovery states are centred; they widen the prediction rather than fitting the 2025 validation mean.

## Run

```matlab
out = run_bridgewest_joint_calibrated_v3_2( ...
    'D:\...\BridgeWest_Stabil_v3_2_JointDiscrepancy', ...
    'D:\...\BridgeWest_Stabil_Shaker_v3_results', ...
    'D:\...\BridgeWest_Joint_results');
```

For a fast check on selected passages:

```matlab
out = run_bridgewest_joint_calibrated_v3_2(packageRoot,resultRoot,outputRoot, ...
    {'2025_run22','2025_run40'});
```

The original `run_bridgewest_stabil_shaker_informed.m` is retained unchanged for comparison.
