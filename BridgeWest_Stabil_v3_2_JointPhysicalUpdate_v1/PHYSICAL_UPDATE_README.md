# Bridge West joint-discrepancy physical FE update v1

## What changed

This branch keeps the established Bridge West joint-discrepancy inference, but replaces the baseline **modal-frequency substitution + deterministic modal-residue eigenvector scaling** with an **equivalent physical update of the FE mass matrix** before inference.

The starting model is the existing Stabil shell-and-beam representation (963 nodes, 864 shells, 188 beams). The raw global modes identified from the shaker spatial shapes are:

| family | raw FE mode | raw frequency (Hz) | shaker spatial MAC |
|---|---:|---:|---:|
| global bending | 3 | 5.2110 | 0.9987 |
| global torsion | 8 | 7.7364 | 0.9961 |

The shaker low-force targets are 6.53534 and 8.76083 Hz. Directly overwriting those eigenvalues reproduces the poles but does not create a physically consistent FE baseline. The new branch updates `M`, re-solves `K phi = omega^2 M phi`, re-identifies the physical modes, and rebuilds all strain/displacement observation operators.

## Physical-update parameterisation

The mass audit of the current reconstructed model gives approximately 24.43 t structural steel, 0.88 t walkway surface mass, and 100.36 t residual nonstructural mass needed to reach the prescribed 6830 kg/m total line mass. Because this residual track/ballast component is not explicitly modelled, the shaker update identifies two *equivalent dynamic* parameters:

1. `eta`: fraction of the residual nonstructural mass participating in the shaker-frequency FE model;
2. `lambda`: fraction of that participating residual mass redistributed transversely, at each longitudinal x-slice, from the tributary deck distribution toward the two girder/deck lines.

The sourced static mass inventory is **not** rewritten. These are equivalent dynamic-mass parameters, not a claim that a measured quantity of ballast physically disappears.

The Python pilot solution is:

- `eta = 0.58501716`;
- `lambda = 0.14232127`;
- effective dynamic total mass = 84.024 t (66.86% of the nominal static total);
- re-solved bending family = FE mode 5 at 6.53530 Hz, MAC = 0.99695;
- re-solved torsion family = FE mode 8 at 8.76080 Hz, MAC = 0.99858;
- force-normalised modal-residue ratios to shaker = 0.969 and 1.109.

Thus the calibrated FE no longer needs deterministic eigenvector/residue scaling to obtain a plausible shaker force-response scale. Modal order is allowed to change: the raw bending family is mode 3, while the corresponding re-solved bending family is mode 5.

## Information policy

The **physical FE update itself uses shaker evidence only**: low-force modal frequencies plus the supplied four-station shaker shapes. It does not read 2022/2025 train-passage strain, static weighing, withheld target strain, or PLS displacement.

To remain compatible with the established joint-discrepancy branch, the runner then retains its existing observation-recovery step: 2022 independent static bogie weighing together with A031/B053 determines one frozen strain-recovery scale. That step is separate from FE physical calibration. The 2025 A031/B053 records then infer passage-specific loads; 2025 static weighing, withheld target strains and PLS displacement are opened only after the requested predictions are frozen.

## Pilot check

The physical update was screened on two 2025 passages without using their validation measurements to fit the FE update.

| event | inferred total (kN) | static (kN) | error | static in 95% total interval | BA peak ratio | BB peak ratio |
|---|---:|---:|---:|---|---:|---:|
| 2025 run 22 | 2059.8 | 2197.4 | -6.3% | yes | 0.926 | 0.925 |
| 2025 run 40 | 2302.2 | 2197.4 | +4.8% | yes | 0.871 | 0.877 |

The mean displacement scale is therefore much closer to the PLS measurements than the original uncorrected physical FE. However, the run40 displacement interval coverage remains poor (about 6.6% BA and 4.8% BB), so the **uncertainty model is not yet uniformly calibrated even though the mean physics improves**.

## MATLAB files added

- `run_bridgewest_joint_physical_v1.m` — complete physical-update + joint-inference runner.
- `matlab/bw_shaker_physical_mass_update.m` — shaker-only mass-matrix calibration and re-solved eigensystem.
- `config/physical_update.json` — parameter bounds, pilot starting point and scientific policy.
- `joint_discrepancy/physical_mass_update_pilot.json` — internal pilot audit.
- `joint_discrepancy/pilot_results/physical_update_run22_run40.csv` — compact two-passage pilot screen.

The original `run_bridgewest_joint_calibrated_v3_2.m` remains unchanged for comparison.

## Run

```matlab
packageRoot = 'D:\BridgeWest_Stabil_v3_2_JointPhysicalUpdate_v1';
resultRoot  = 'D:\BridgeWest_Stabil_Shaker_v3_results';
outputRoot  = 'D:\BridgeWest_JointPhysical_results';

out = run_bridgewest_joint_physical_v1(packageRoot,resultRoot,outputRoot);
```

For a two-passage check:

```matlab
out = run_bridgewest_joint_physical_v1( ...
    packageRoot,resultRoot,outputRoot, ...
    {'2025_run22','2025_run40'});
```

Then plot the frozen results with the same joint plotter:

```matlab
figs = plot_bridgewest_joint_results( ...
    packageRoot,resultRoot,outputRoot, ...
    'Event','all','Visible','off');
```

## Important interpretation

This is a **stronger physical baseline than direct frequency substitution** because the FE matrix is changed first and the eigensystem and response operators are re-solved. It is still an *equivalent* physical update: two shaker families cannot uniquely identify literal ballast participation, support stiffness, every section property and all connection parameters. The calibrated parameters should therefore be reported as effective dynamic mass-distribution parameters of the reconstructed FE model.
