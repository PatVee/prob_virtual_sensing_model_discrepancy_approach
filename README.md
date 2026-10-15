# Probabilistic virtual sensing — quick start

MATLAB scripts for **Bridge West** (shaker-informed physical FE update and virtual
sensing) and **KW51** (blueprint FE ensemble and modal-only discrepancy).
MATLAB with Signal Processing Toolbox is required; bundled FE/Stabil code is included.

## Bridge West

Main analysis: `BridgeWest_Stabil_v3_2_JointPhysicalUpdate_v1/run_bridgewest_joint_physical_v1.m`.
An existing `BridgeWest_Stabil_Shaker_v3_results` folder with prepared 2022/2025
passages is required (not included in this archive).

```matlab
bwRoot = 'C:\path\BridgeWest_Stabil_v3_2_JointPhysicalUpdate_v1';
v3Root = 'C:\path\BridgeWest_Stabil_Shaker_v3_results';
outRoot = fullfile(bwRoot,'joint_physical_results');
addpath(bwRoot);
out = run_bridgewest_joint_physical_v1(bwRoot,v3Root,outRoot);
plot_bridgewest_joint_results(bwRoot,v3Root,outRoot,'Event','all','Visible','off');
```

The Bridge West `main.m` only builds/plots the Lorenzen FE model; it is **not**
the full prediction runner. If needed, first generate v3 results using
`production/matlab/run_all_bridgewest_stabil_v3.m` with your raw-data paths.

## KW51

**Recommended order when documenting both methods:**

1. **Original blueprint ensemble (optional):** put the original monitoring ZIPs
   in `KW51_Probabilistic_Ensemble_MATLAB_BLUEPRINT_R2025a_FIX/KW51_Probabilistic_Ensemble_MATLAB/data/`,
   then run that folder's `main.m`.
2. **Modal-only approach:** edit `zipRoot` in
   `KW51_modal_only/main_kw51_modal_only.m` to point to the
   `KW51_Probabilistic_Ensemble_MATLAB_BLUEPRINT_R2025a_FIX` folder; then run
   the following **from `KW51_modal_only/`**:

   ```matlab
   main_kw51_modal_only
   ```

The required raw inputs are `traindata_201912_part1.zip` and
`traindata_201912_part2.zip` (not bundled). The verified modal bank **is** bundled.
**The original `main.m` need not finish—or even run—before modal-only**:
`main_kw51_modal_only.m` independently prepares raw inputs and uses the fixed,
unscaled blueprint FE with measured frequency/damping variability, not the earlier
sampled-stiffness ensemble.

Modal results: `KW51_Probabilistic_Ensemble_MATLAB_BLUEPRINT_R2025a_FIX/KW51_modal_only_results/modal_YYYYMMDD_HHMMSS/`.
Replot saved results without rerunning inference:

```matlab
plot_kw51_modal_only(out,'Visible','on');
% or: plot_kw51_modal_only('C:\path\modal_YYYYMMDD_HHMMSS','Visible','on');
```

