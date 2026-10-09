> Historical v3/v3.1 document retained for provenance. For the current update, read `V3_2_README.md` and `audit/MASS_AND_LOAD_AUDIT.md`. No mass correction is activated in v3.2.

# v3.1 implementation summary

The production workflow now implements the requested formulation directly:

1. The detailed Stabil model and 84-state shaker-informed discrepancy ensemble are frozen.
2. The configured event set contains 12 Bridge West passages from 2022 and 24 from 2025.
3. For each passage, only A031 and B053 structural strains are read for inference. Speed and bridge-entry time use the independent upstream rail strain channels.
4. A bounded continuous posterior is inferred for the eight bogie resultants.
5. The shaker-informed discrepancy is conditioned on A031/B053 and propagated to A058, A070, B031 and B081.
6. The prediction is saved before any target measurement is opened.
7. A separate scoring stage reopens the raw passage and reads only the four target channels using the already-frozen speed and entry time.
8. All four target predictions are retained for every successfully predicted passage. Non-`ok` target observations remain documented but do not enter the primary aggregate.
9. No second target-derived discrepancy model is fitted.

A selective HDF5 reader was added so inference does not load the entire bridge signal matrix merely to select A031/B053. The kinematic routine likewise reads only the upstream rail timing channels.

`production/verify_v3_package.py` performs static checks on event counts, channel configuration, selective reading, target-stage separation and obsolete discrepancy wording.


## Displacement plot implementation

The prediction MAT file now saves the exact 41-point dataset-coordinate vector from `bw_lorenzen_observation_maps` together with BA/BB displacement means and nominal 95% intervals. `score_bridgewest_stabil_displacement.m` opens 2025 PLS records only after prediction, computes independent measured/predicted own-midspan-peak spatial profiles, and writes displacement score/profile tables. `plot_bridgewest_stabil_displacement.m` creates (i) the selected-event BA/BB spatial displacement figure with PLS overlay when available and (ii) a 2025 campaign summary of peak-ratio and spatial NRMSE.
