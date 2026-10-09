> Historical v3/v3.1 document retained for provenance. For the current update, read `V3_2_README.md` and `audit/MASS_AND_LOAD_AUDIT.md`. No mass correction is activated in v3.2.

# Frozen 36-passage reference-only protocol

## Research calculation

The Bridge West Stabil model and the 84-state shaker-informed discrepancy ensemble are frozen before train-passage inference.

For every selected train passage, the structural measurements permitted during inference are exactly:

- `A031` (`str_x_BA_uf_031`)
- `B053` (`str_x_BB_uf_053`)

The default event set contains 36 passages:

- 12 passages from 2022: 02, 04, 06, 08, 10, 12, 14, 16, 18, 20, 24, 26;
- 24 passages from 2025: 02, 04, 06, 08, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50.

For each passage, A031/B053 are used to infer the bounded eight-bogie load posterior and condition the shaker-informed Stabil discrepancy. The conditioned structural response is propagated to four withheld strain locations:

- A058
- A070
- B031
- B081

The target measurements do not enter the prediction stage. They are first read in the separate assessment function after the corresponding prediction MAT file already exists.

## Strict software separation

`run_bridgewest_stabil_shaker_informed` preprocesses only A031 and B053. Its saved `permitted_inputs/<event>.mat` contains two strain columns only. It does not read or preprocess A058, A070, B031 or B081.

`score_bridgewest_stabil_results` is a separate post-prediction stage. It reopens the raw passage, uses the already-frozen speed and entry time from prediction, and extracts only the four protected target channels for assessment. Target data cannot change inferred loads, the Stabil model, shaker-state weights or predictive intervals.

All four targets are predicted for all successfully inferred passages. A target measurement with a non-`ok` quality flag remains in `target_assessment_all.csv` when numerically available, but it is excluded from the primary aggregate table. This retains the 2022 A070 predictions without treating reduced-quality observations as primary validation evidence.

No second target-derived discrepancy model is fitted in this protocol.


## Displacement prediction and validation

For every successfully inferred passage, the same A031/B053-conditioned posterior is propagated to 41 vertical-displacement stations on BA and 41 on BB. No displacement measurement is used to infer the bogie loads or to condition the shaker-informed discrepancy.

For 2025 passages with PLS data, `score_bridgewest_stabil_displacement` is run only after the prediction MAT files already exist. It extracts an independent measured own-midspan-peak spatial profile for each girder and compares it with the predicted own-midspan-peak profile. The PLS profiles are assessment-only and cannot alter the FE model, load posterior, shaker-state weights or displacement intervals.

The displacement plots therefore show an independent cross-response check of the strain-conditioned Stabil model rather than an additional calibration step.
