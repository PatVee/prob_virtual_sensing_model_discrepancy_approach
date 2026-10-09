# KW51: modal-discrepancy-only workflow

Use this folder as a separate add-on beside the original
`KW51_Probabilistic_Ensemble_MATLAB` package. Existing scripts and results are not
edited. There is one prediction method and one output folder per run.

## Run

1. Extract this ZIP. Open `KW51_modal_only/main_kw51_modal_only.m`.
2. Keep or adjust `zipRoot` to your existing blueprint-package location. The
   script already contains the path supplied in this conversation.
3. Confirm the two ZIP filenames in `originalRoot/data`. Change them to
   `_copy.zip` only if that is how your files are named.
4. Run `main_kw51_modal_only`.

Do not run the previous comparison main. This main calls `run_kw51_modal_only`
and uses the included verified modal bank; a native FE rebuild is not required.
The original package is a dependency for input preparation, FE load maps,
kinematic fitting, Gaussian conditioning and target scoring.

Direct call:

```matlab
addpath('C:\path\KW51_modal_only');
out = run_kw51_modal_only(originalRoot,passageZIPs,outputRoot);
```

To start with the two previously selected passages:

```matlab
out = run_kw51_modal_only(originalRoot,passageZIPs,outputRoot, ...
    'Events',{'traindata_20191215_081100','traindata_20191217_080627'});
```

Each execution creates `KW51_modal_only_results/modal_YYYYMMDD_HHMMSS`.
The runner saves `run_output.mat` before plotting, so plotting can be repeated
without repeating inference:

```matlab
plot_kw51_modal_only(out,'Visible','on');
% Alternatively:
plot_kw51_modal_only('C:\path\modal_YYYYMMDD_HHMMSS','Visible','on');
```

Separate discrepancy diagnostics, including optional plots:

```matlab
quantify_kw51_modal_discrepancy(out,'Plot',true,'Visible','on');
```

The output struct now has `run_dir`, not separate method folders. Do not call
`run_kw51_from_archives`, `run_kw51_modal_comparison`, or the old comparison plot.
Independent kinematic-truth plotting is not automatically requested: it requires
an independently supplied truth file and is separate from strain validation.

## Model now executed

For every candidate axle-load vector, the nominal response is the same unscaled
blueprint quasi-static FE operator applied to those loads. There is no sampled
stiffness state and no division by 0.8, 0.9 or 1.0 in the nominal ensemble.
The fixed stiffness scale is exactly one.

Modal discrepancy uses contemporaneous joint rows of the measured frequency and
damping records. For matched mode r:

- baseline frequency = unscaled FE frequency;
- perturbed frequency = baseline frequency × measured frequency / measured median;
- baseline damping = measured median damping ratio;
- perturbed damping = sampled measured damping ratio.

Each paired perturbation uses the same axle loads and fitted kinematics in the
perturbed and baseline dynamic responses. Their difference supplies the
cross-output discrepancy covariance. Reference strains condition the candidate
loads and discrepancy using the original Gaussian-mixture calculation.

The six tracked modes 1–6 map to FE modes 1,2,3,5,4,6. The supplied pre-December
window supplies 916 complete joint samples. Percent damping was converted to a
ratio during preparation. No mode-selection threshold was changed for this update.
Every passage is checked to occur after the last modal sample.

Removed from the executed model:

- discrete nominal stiffness states and the 0.9 paired baseline;
- prescribed random stiffness/compliance and strain-recovery factors;
- prescribed support-force and kinematic perturbations;
- support-force covariance weighting in the kinematic fit;
- all method-comparison runs and plots.

Retained assumptions include the load prior, fitted kinematics and their fitting
regularisation, measurement/output noise, covariance regularisation, likelihood
stride/power, fixed mode shapes, the selected mode subset and empirical window.
Kinematics are fitted once from permitted reference/rail data and are not sampled.
Measured modal scatter includes environmental variation and identification error.

**This is modal-informed discrepancy around a fixed FE baseline.** It does not
calibrate physical FE parameters to the measured mean frequencies, and it does
not replace the nominal quasi-static prediction with a fully dynamic prediction.
It does not claim exact Bridge West implementation parity without that source.

## Outputs and interpretation

- `predictions/`, `draws/`, `decomposition/`: one method only.
- `validation_scores.csv` and `validation_summary.csv`: withheld-target scores.
- `figures/`: measured histories, modal-only means and nominal 95% intervals,
  exported as PDF and PNG.
- `discrepancy_diagnostics/`: posterior-weighted FE contribution, inferred
  correction and exact pointwise variance decomposition.
- `validation_skipped.csv`: reasons why some measured targets were unavailable.

The FE contribution and correction use the same posterior load weights. Their
comparison is a decomposition, not a separately fitted competing model. The
inferred correction is not independently identified true model discrepancy.
The FE variance now represents load variation only. Covariance terms may be
negative and are not independent percentage shares. Raw variance components
are unfiltered; history diagnostics use the same filtering/baseline convention
as validation.

Predictions, draws, decomposition, input files, FE/modal assets and inference code
are frozen before scoring protected targets. The add-on's plotting and diagnostic
functions are excluded from its inference-code hash list so display fixes can be
made without altering saved inference. Do not modify frozen numerical functions
for an existing run; start a new run after such changes.

## Executed checks and a material limitation

The complete numerical runner, including validation and CSV diagnostics, ran in
GNU Octave 8.4.0 on both selected passages with 900 candidates, 48 paired draws
and 256 posterior draws. Modal resolution and decomposition checks passed, as did
unit checks for the unscaled modal centre and missing validation channels.
All original package files remain unchanged.

**Both test passages have posterior ESS approximately one, and target coverage is
low.** Removing the broader perturbations produces narrow intervals that do not
cover the observed errors well. This is an observed scientific limitation, not a
runtime error. These results do not justify a claim of improved uncertainty
calibration or replacing physical FE calibration with measured modal scatter.
No settings were tuned using the withheld targets to repair coverage.

See `execution_checks/` for actual test scores and numerical checks. These are
two-passage development checks, not a completed full-campaign study. Native MATLAB
plotting has not been executed here; data-shape/missing-channel checks were run.
MATLAB R2025a and Signal Processing Toolbox remain the intended plotting runtime.

## Optional checks

```matlab
test_kw51_modal_only(originalRoot);
```

`'RebuildModalBank',true` rebuilds the FE bank with the previously audited scaled
eigensolver and static replay checks. The default cached-bank route verifies its
FE source hashes and exact static response arrays against the original package.
`'PointsPerCycle',80` increases integration resolution. The first passage is
always checked against doubled resolution; `'CheckEveryEvent',true` checks all.
