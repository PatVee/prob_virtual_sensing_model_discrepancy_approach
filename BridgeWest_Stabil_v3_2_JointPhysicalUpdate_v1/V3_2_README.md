# Bridge West v3.2: mass/load audit and restored strain plotting

## What this update does

The strain plotter is repaired. The v3.1 wrapper plotted displacement only; the same command now plots the four withheld strain targets for the selected passage. Selection supports 2022, 2025, both campaigns, and groups of runs. Displacement remains accessible through an explicit option.

The audit independently verifies the static train weighing and the mass bookkeeping in the supplied native results. **It does not justify a new physical bridge mass. No mass, stiffness, load-posterior or displacement multiplier has been applied.** Read `audit/MASS_AND_LOAD_AUDIT.md` for the findings.

The original numerical inference remains A031/B053-only, with A058/A070/B031/B081 and PLS reserved for assessment. The optional known-static-load diagnostic is a separately named forward comparison, not a hidden change to this protocol.

## Install

This complete package contains the original v3 files, the v3.1 displacement additions, and the v3.2 plotting/audit update. Alternatively, the small patch ZIP can be overlaid onto the root of an existing complete v3/v3.1 package. Preserve its `production`, `data`, `tests` and `audit` folder structure.

The older v3.1 download is an overlay, not a complete stand-alone FE package. Use this complete package if the original full v3 installation is no longer available.

In MATLAB, set the path to the chosen installation:

```matlab
packageRoot = 'D:\BridgeWest_Stabil_v3_2_MassAudit_Plots';
addpath(fullfile(packageRoot,'production','matlab'),'-begin');
clear plot_bridgewest_stabil_results bw_select_saved_events bw_stabil_strain_plot_data
rehash
which plot_bridgewest_stabil_results -all
```

The first path shown must be the updated file. Older same-name functions elsewhere on the MATLAB path can otherwise hide the repair.

**Do not rerun the FE solver merely to plot existing results.** Set `outputRoot` to the existing result folder containing `predictions`, `validation_targets`, and `tables`.

## Original single-event command

```matlab
plot_bridgewest_stabil_results(outputRoot, ...
    'DataRoot',dataRoot,'Event','2025_run22');

plot_bridgewest_stabil_results(outputRoot, ...
    'DataRoot',dataRoot,'Event','2022_run20');
```

Cached `validation_targets/<event>.mat` files are used first. Thus `DataRoot` is optional for the supplied completed result archive. When a cache is missing, a supplied `DataRoot` allows assessment-only target extraction using the saved passage speed, entry time and prediction grid. Without either source, the four predictions are still plotted with measurements labelled unavailable.

## Plot complete campaigns or selected groups

```matlab
% All 36 saved passages, no open windows left after saving:
plot_bridgewest_stabil_results(outputRoot,'Event','all', ...
    'Visible','off','CloseAfterSave',true);

% All saved 2022 passages:
plot_bridgewest_stabil_results(outputRoot,'Event','2022', ...
    'Visible','off','CloseAfterSave',true);

% All saved 2025 passages:
plot_bridgewest_stabil_results(outputRoot,'Event','2025', ...
    'Visible','off','CloseAfterSave',true);

% Three particular 2025 passages:
plot_bridgewest_stabil_results(outputRoot,'Campaign',2025,'Runs',[20 22 40]);

% A mixed-campaign selection:
plot_bridgewest_stabil_results(outputRoot, ...
    'Event',{'2022_run20','2025_run22'});
```

The function selects **saved predictions**, not arbitrary raw runs. The configured 2022 list contains runs 2,4,6,8,10,12,14,16,18,20,24,26: there is no saved 2022 run22. Missing requested predictions cause an explicit error showing the available events; the plotter never substitutes another campaign or passage.

## Outputs and quality handling

Each event generates a four-panel A058/A070/B031/B081 strain figure and saves PNG, vector PDF and editable FIG files under `outputRoot/figures`:

    bw_stabil_strain_2025_run22.png / .pdf / .fig
    bw_stabil_strain_2022_run20.png / .pdf / .fig

Titles state the campaign/run and the same range-normalised RMSE, correlation and coverage definitions as the saved scores. Nominal 95% bands are not retuned. Non-ok target measurements are not scored; the prediction remains visible. `ShowReducedQuality=true` permits an explicitly labelled dotted reduced-quality overlay, still excluded from primary scores. No toolbar should appear in exported figures.

Optional settings: `Save`, `Visible`, `OutputDir`, `FontSize`, `DPI`, `CloseAfterSave`, `ShowReducedQuality`.

## Retain displacement plots

```matlab
% Displacement only:
plot_bridgewest_stabil_results(outputRoot,'Event','2025_run22', ...
    'Figures','displacement');

% Both strain and displacement:
plot_bridgewest_stabil_results(outputRoot,'Event','2025_run22', ...
    'Figures',{'strain','displacement'});
```

Displacement uses saved prediction fields and assessment-only PLS snapshots. It does not fit a gain to PLS or silently rerun PLS scoring. To regenerate missing PLS assessment files, explicitly use the existing `score_bridgewest_stabil_displacement(dataRoot,outputRoot)`.

The existing displacement scores and their missing-profile handling have not been changed by this update. The documented 39-profile result is not silently relabelled as a larger assessment.

## Run the mass/load audit locally

```matlab
A = audit_bridgewest_mass_and_load(outputRoot);
```

This reads the saved FE matrix and posterior loads, checks rigid-translation mass and component bookkeeping, and compares inferred totals with the original static weighing. Outputs go to `outputRoot/audit_mass_load`. It does not open target strains or PLS and does not modify predictions.

For an exact fixed-FE load comparison, rather than a final-curve multiplier:

```matlab
run_bridgewest_known_load_diagnostic(outputRoot,'Event','2025_run22');
% Or: run_bridgewest_known_load_diagnostic(outputRoot,'Event','all');
```

This optional function regenerates the nominal unit-bogie basis and combines it with the saved inferred loads and the actual weighed bogie pattern separately. Both strain and displacement are saved under `known_load_diagnostic`. It does not recompute reference conditioning or uncertainty. It is not a corrected production prediction and has not been numerically executed here in MATLAB.

## Tests and status

```matlab
addpath(fullfile(packageRoot,'tests'));
test_bridgewest_stabil_v3_2(outputRoot);
```

This native-MATLAB test checks selections, channel mappings, quality flags, 132 primary scores and the two weighing totals without requiring graphics. It expects the supplied complete 36-passage result archive.

Executed during preparation: Python reading of native v7/v7.3 MATLAB results, mass integration, independent weighing totals, score recomputation for 144 histories, and before/after hashes of all 36 predictions. Maximum score differences are below 2e-15. **MATLAB execution, graphics export and a new FE campaign run were not available/performed here.** The executable audit is `tests/verify_saved_results.py`; `audit/verification.json` records exactly what was checked.

Legacy model reports remain for traceability. Their inherited numerical assumptions must not be treated as new independent mass evidence.
