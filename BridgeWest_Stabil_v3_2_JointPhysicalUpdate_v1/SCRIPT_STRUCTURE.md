> Historical v3/v3.1 document retained for provenance. For the current update, read `V3_2_README.md` and `audit/MASS_AND_LOAD_AUDIT.md`. No mass correction is activated in v3.2.

# Script structure

## Detailed Stabil model

| Stage | MATLAB file | Use |
|---|---|---|
| FE entry point | `main.m` | Runs the literature-based detailed model and modal plots |
| Configuration | `matlab/run_bridgewest_lorenzen.m` | Reads model JSON, solves and saves the detailed FE state |
| Mesh | `matlab/bw_lorenzen_mesh.m` | Builds shell/beam nodes, connectivity and bearing constraints |
| Assembly | `matlab/bw_lorenzen_assemble.m` | Assembles steel, added mass and Garather-Bach upper-flange restraint |
| Modal solve | `matlab/bw_lorenzen_modes.m` | Solves and mass-normalises modes |
| Observation maps | `matlab/bw_lorenzen_observation_maps.m` | Exact shell surface-strain maps and BA/BB displacement maps |
| Shaker ensemble | `shaker_informed/matlab/bw_shaker_ensemble.m` | Loads the frozen 84 weighted modal states |
| Passage basis | `shaker_informed/matlab/bw_stabil_basis_shaker.m` | Generates unit-bogie response bases for a passage speed/state |
| Paired discrepancy | `shaker_informed/matlab/bw_stabil_paired_discrepancy.m` | Generates nominal-minus-shaker-state response differences |

## 36-passage production protocol

| Stage | MATLAB file | Use |
|---|---|---|
| Run all | `production/matlab/run_all_bridgewest_stabil_v3.m` | Executes prediction, audit, assessment and plots |
| Reference-only prediction | `production/matlab/run_bridgewest_stabil_shaker_informed.m` | Reads only A031/B053 structural histories for each passage and predicts four targets |
| Load posterior | `production/matlab/bw_continuous_load_posterior.m` | Infers the bounded eight-bogie posterior from A031/B053 |
| Shaker-conditioned prediction | `production/matlab/bw_predict_continuous_shaker.m` | Conditions shaker discrepancy on A031/B053 and propagates uncertainty |
| Strict inference preprocessing | `production/matlab/raw_functions/bw_prepare_raw_event.m` | Explicitly requires a permitted channel list; production passes A031/B053 only |
| Post-prediction target extraction | `production/matlab/raw_functions/bw_prepare_validation_targets.m` | Opens A058/A070/B031/B081 only after prediction |
| Strain assessment | `production/matlab/score_bridgewest_stabil_results.m` | Calculates target errors/coverage; never feeds results back to inference |
| Displacement assessment | `production/matlab/score_bridgewest_stabil_displacement.m` | Opens 2025 PLS only after prediction and scores BA/BB own-peak spatial profiles |
| PLS snapshot reader | `production/matlab/raw_functions/bw_pls_raw_snapshot.m` | Assessment-only extraction of measured BA/BB own-midspan-peak profiles |
| Firewall audit | `production/matlab/audit_bridgewest_reference_only_protocol.m` | Verifies two-channel permitted inputs and no target-use flags across 36 events |
| Plot wrapper | `production/matlab/plot_bridgewest_stabil_results.m` | Calls the displacement plotting stage from saved predictions/assessment |
| Displacement plotting | `production/matlab/plot_bridgewest_stabil_displacement.m` | Event BA/BB profiles plus 2025 amplitude-transfer and spatial-NRMSE summary |

The production protocol contains no second target-derived discrepancy calibration.
