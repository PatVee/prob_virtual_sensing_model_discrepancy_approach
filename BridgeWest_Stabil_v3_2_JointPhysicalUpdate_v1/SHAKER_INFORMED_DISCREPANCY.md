> Historical v3/v3.1 document retained for provenance. For the current update, read `V3_2_README.md` and `audit/MASS_AND_LOAD_AUDIT.md`. No mass correction is activated in v3.2.

# Bridge West Stabil model: shaker-informed discrepancy ensemble

## Status

This revision connects the provisional Lorenzen-inspired shell-and-beam Stabil
model to the September 2024 / January 2025 Bridge West shaker tests without
reusing the earlier arbitrary independent FE perturbation standard deviations.
The purpose is to construct the **model-discrepancy ensemble**, not to refit the
bridge to withheld train strains or PLS displacement.

The physical reconstruction is still provisional because the original Garather
Bach SOFiSTiK `.cdb`/`.dat` model and dimensioned plate/stiffener drawings are
not present. The supplied `sofi-service-4-zekiss` archive is a generic SOFiSTiK
CDB interface. It is useful as an extraction route if the original CDB becomes
available, but it does not itself contain the bridge model.

## Lorenzen spring clarification

Lorenzen's **Table 5.2** belongs to the HUMVIB footbridge in Chapter 5. Those
support/connection spring values are therefore **not transferred to Garather
Bach**. In Garather Bach Section 7.2.3, the relevant restraint is the spring
representation at the ends of the upper flanges, with a modulus of subgrade
reaction of **700,000 kN/m3**. That value remains in the reconstructed bridge
model. Its exact footprint and direction remain reconstruction assumptions until
a native model or drawings are available.

## What the shaker data identify

The archive contains 113 forced-vibration records. The modal-identification
audit produces 72 accepted single-mode fits in the BA/BB measurement family.
The detailed Stabil rebuild has a first global bending candidate at FE index 3
and a first global torsional candidate at FE index 8. At the four available
shaker stations their spatial MAC values are approximately 0.9987 and 0.9961,
respectively. The mode shapes are therefore retained in the discrepancy ROM.

The raw reconstructed FE frequencies (about 5.211 and 7.736 Hz) are not treated
as calibrated frequencies. Instead the central low-force shaker state is used
as the nominal modal state:

- first global bending: 6.53534 Hz;
- first global torsion: 8.76083 Hz;
- bending damping ratio: 0.01634;
- torsional damping ratio: 0.01377.

These are the medians of all accepted 5 m/s2 fits for the respective mode. The
quantity being updated is the **modal ROM**, not Young's modulus and not a unique
set of physical section stiffnesses.

## Joint ensemble construction

The two campaigns have different experimental layouts. In 2024, bending and
torsional modes were identified in the same runs, giving 20 genuine accepted
joint states. These same-run pairs preserve the observed cross-mode dependence.

In 2025 the two modal families were excited in separate run blocks. The code
therefore does **not** rank-pair them. Within each shaker acceleration level,
the 2025 joint uncertainty is represented by the product of the two empirical
mode marginals. This is conditional independence, not fabricated same-run
correlation.

The production mixture contains 84 modal states:

- 20 genuine 2024 same-run states;
- 64 2025 within-level product states.

The state weights are explicit. Each campaign has total weight 0.5; within each
campaign the 5, 7.5 and 10 m/s2 excitation levels each have one-third of the
campaign weight; states within a level share that level weight. These weights
are a transparent uncertainty-mixture design, not an estimated long-term
probability distribution of bridge states.

Across the weighted ensemble the observed ranges are approximately:

| quantity | range |
|---|---:|
| bending frequency | 6.3335--6.5402 Hz |
| torsional frequency | 8.5217--8.8691 Hz |
| bending damping ratio | 0.01552--0.02124 |
| torsional damping ratio | 0.01318--0.02441 |
| bending eigenvalue scale relative to nominal low-force state | 0.9392--1.0015 |
| torsional eigenvalue scale relative to nominal low-force state | 0.9462--1.0249 |

The 2025 torsional fits show a clear campaign shift relative to 2024. This is
why the final mixture includes both campaigns instead of treating the 2024
same-run covariance as the complete uncertainty description.

## Modal-state propagation through the Stabil model

For state `m`, the detailed FE mode shapes, mass distribution, geometry and
observation maps are fixed. Only the two experimentally identified modal
frequencies and the Rayleigh-damping anchors are changed:

```
omega_1,m = 2*pi*f_1,m
omega_2,m = 2*pi*f_2,m
```

with Rayleigh coefficients selected so that the measured damping ratios are
reproduced at those two modes. Other retained modes use the corresponding
Rayleigh curve. This is intentionally a modal update: the shaker tests do not
uniquely identify separate physical values for Young's modulus, plate thickness,
ballast mass, cross-member stiffness and spring footprint.

The dynamic runner uses mass-normalised Stabil modes up to 25 Hz by default.
A 16-axle / 8-bogie moving-load basis is applied as vertical axle resultants on
the deck track centreline. This avoids inventing unverified rail coordinates.
The existing axle offsets are retained. Shell-surface strain is recovered by
the Stabil stress operator; the old fitted beam curvature-to-strain recovery
depth is not used.

## Paired discrepancy definition

The implementation first forms an eight-bogie **unit-load discrepancy basis** for each shaker state:

```
DeltaB_m(t) = B(thetam) - B(theta0)
```

where `theta0` is the nominal low-force shaker state and `thetam` is one of the weighted shaker states. For any passage load vector `p`,

```
delta_m(t | p) = DeltaB_m(t) p.
```

This is stronger than assigning a different random load to every discrepancy state: the structural/modal discrepancy is stored independently of load uncertainty and can be projected with a passage-specific candidate, posterior mean, or explicit load quadrature later.

The response blocks are

1. seven structural strain channels in microstrain;
2. 41 BA vertical displacement stations in mm;
3. 41 BB vertical displacement stations in mm.

The reference channels are A031 and B053, corresponding to strain indices 1
and 4 in the new observation map. The weighted discrepancy moments are then
computed from the 84 state responses. Reference-target cross-covariance supplies
the conditional discrepancy correction used by the virtual-strain predictor.

## Important change from the old ensemble

Do **not** simultaneously apply the previous independent perturbation standard
deviations for rigidity, mass, transverse coupling, end restraint, damping,
common stiffness-mass scale and strain recovery. Doing so would count much of
the same structural uncertainty twice.

For the new Stabil branch, the default separation is:

- uncertain train loading -> candidate/load library;
- shaker-observed structural dynamics -> weighted discrepancy ensemble;
- measurement noise -> separate observation term;
- local strain-transfer error not explained by the shell model -> diagnosed
  from held-out field residuals, not silently absorbed by changing FE parameters.

The old global stiffness candidate grid should therefore not be interpreted as
an additional physical uncertainty prior in this branch. If retained for a
specific sensitivity study, it must be reported as a distinct second layer.

## MATLAB workflow

Run the base Stabil model first and then the discrepancy functions:

```matlab
addpath('matlab');
addpath('shaker_informed/matlab');
fe = run_bridgewest_lorenzen();

root = fullfile(pwd,'shaker_informed');
c = jsondecode(fileread(fullfile(root,'discrepancy_config.json')));
ens = bw_shaker_ensemble(root);
```

For a passage with known speed `v`, observation time vector `tobs`, and an
`M x 8` set of perturbation bogie loads in N:

```matlab
out = bw_stabil_paired_discrepancy(fe,v,tobs,[],c,+1);
% Project the structural discrepancy basis using a passage load vector pN:
delta = bw_project_discrepancy_basis(out.delta_basis,pN);
refs = c.reference_indices_matlab;
disc = bw_simulation_discrepancy_weighted(delta,out.weights,refs,c);
```

For the nominal load-candidate library:

```matlab
[nom,B0] = bw_stabil_nominal_candidates(fe,v,tobs,candidateLoadsN,c,+1);
```

The new predictor `bw_predict_shaker.m` accepts these nominal candidates, the
paired discrepancy states and their weights. Target strains must remain absent
from conditioning.

`run_shaker_informed_discrepancy_demo.m` provides an end-to-end mechanics demo
using reproducible synthetic bogie loads only. It does not calibrate to withheld
train strain or PLS data.

## SOFiSTiK CDB upgrade path

The attached service package exposes the SOFiSTiK 2020 CDB records needed to
upgrade the reconstruction. `shaker_informed/sofi_cdb_adapter` contains an
adapter and instructions. If the original Garather Bach CDB becomes available,
extract its nodes, QUAD thickness/material records, beam section references and
elastic spring records. Replace provisional `model.json` properties with those
source values, rerun the modal mapping, and regenerate the paired discrepancy
responses. The statistical construction of the shaker ensemble does not need
to change.

## Numerical pilot check completed in this package

A Python counterpart was used only to verify the mechanics of the modal-state
propagation because MATLAB/Stabil execution is not available in the delivery
environment. At 160 km/h with a common 270 kN load assigned to each of the eight
bogie bases, the weighted shaker ensemble produced finite, non-zero midspan
vertical-displacement discrepancies. The nominal peak displacement magnitudes
were approximately 4.18 mm (BA) and 4.31 mm (BB). The weighted mean peak
absolute discrepancy was approximately 0.264 mm (BA) and 0.235 mm (BB), and the
maximum pointwise weighted standard deviation was approximately 0.133 mm on
both girders. These are mechanics checks, not field-validation scores.

The pilot outputs are in:

- `shaker_informed/results/python_dynamic_pilot_summary.json`;
- `shaker_informed/results/python_dynamic_pilot_discrepancy.csv`;
- `shaker_informed/results/python_dynamic_pilot_discrepancy.png`;
- `shaker_informed/results/shaker_frequency_state_cloud.png`.

The final strain discrepancy should still be generated in MATLAB using the
Stabil `se_shell4` observation map, because that is the production strain
recovery path in this model.

## Displacement plotting in v3.1

The same shaker-conditioned posterior used for the withheld strains is propagated through the detailed 41-point BA and 41-point BB displacement observation maps. The prediction files save the exact dataset-coordinate vector. PLS measurements, where available in 2025, are opened only after prediction for an independent own-midspan-peak spatial-profile comparison. They do not condition the train-load posterior or the shaker discrepancy.
