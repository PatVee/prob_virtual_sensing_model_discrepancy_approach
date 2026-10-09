# Bridge West: independent mass and train-load audit

## Decision

**No independently supported new bridge mass was established from the available sources. No 0.61 mass/stiffness scale, displacement-only multiplier, or other amplitude correction has been applied.**

The train weighing tables are usable independent information. The bridge's current aggregate mass is a numerical input, not a verified mass inventory. The defensible immediate update is therefore to preserve the unmodified reference-only prediction as the baseline, supply the verified weighing data and a separate known-static-load diagnostic, and repair plotting. A physically supported bridge-mass correction remains conditional on additional evidence.

This audit uses the uploaded native `BridgeWest_Stabil_Shaker_v3_results` archive. It does not substitute the separate effective Timoshenko calibration branch for the current detailed shell-and-beam model.

## 1. What the original sources establish

Lorenzen (2023), Chapter 7, describes the Garather Bach bridge. Figure 7.2 (printed page 103; PDF page 131) establishes the 16.40 m bearing span and the 1 m extensions at each end. This is a sensor/layout drawing, not a complete plate-thickness, cross-section, ballast or sleeper schedule.

Section 7.2.3 (printed page 111; PDF page 139) says that the original SOFiSTiK model was built from the design plans. The main load-bearing structure uses shell elements, the walkway and handrail use beam elements, and ballast is represented as additional mass. Springs at the upper-flange ends represent partial restraint from the track and ballast. Their reported subgrade-reaction modulus is 700,000 kN/m^3. The passage does not supply the mass quantities needed to reconstruct ballast, rails, sleepers and ancillary dead load independently.

The explicit cross-reference to Table 5.1 concerns steel properties. Table 5.2 describes the HUMVIB footbridge's springs; it is not a Bridge West spring inventory.

The supplied `sofi-service-4-zekiss-main(1).zip` is a generic interface. It includes `dat/example_frame.dat` and `testing/test.dat`, which contain the same four-node example frame. It does not contain a Garather Bach CDB or bridge-specific input model from which its mass could be integrated. The presence of example `.dat` files must not be confused with possession of the bridge model.

The recovered measurement metadata contain the original weighing protocols and sensor positions. They do not provide the missing structural mass bill.

## 2. What the saved model actually contains

The model configuration imposes

    nominal line mass = 6,830 kg/m
    model length      = 18.40 m
    imposed total     = 125,672 kg.

The original `parameter_provenance.csv` labels the line mass with `Broetzmann2025 p152`. That citation could not be verified from a matching original source. Its value is retained as an inherited nominal parameter, not promoted to an independently measured mass. The original provenance table is preserved in this audit folder; the distributed table now explicitly marks the line-mass source as unverified.

Recalculation from the uploaded native mass matrix gives:

| Quantity | Mass (kg) | Share of total |
|---|---:|---:|
| Steel from the numerical elements | 24,427.406969 | 19.4374% |
| Assigned walkway surface mass | 883.200000 | 0.7028% |
| Residual added mass | 100,361.393031 | 79.8598% |
| **Total** | **125,672.000000** | **100%** |

These are **computed model quantities**, not three separately measured component masses. Steel and walkway estimates still depend on provisional dimensions/densities.

The decisive assembly operation in `matlab/bw_lorenzen_assemble.m` is:

```matlab
residual = c.mass_per_length_total_kg_m*c.length_m ...
           - steelMass - sum(walkMass);
added = walkMass + residual*deckArea/sum(deckArea);
```

Thus the residual added mass closes a prescribed total. It is not a ballast volume multiplied by a documented ballast density. Increasing or decreasing assumed steel density alone reallocates the residual, while the imposed aggregate mass remains fixed.

The independent numerical check uses rigid translation:

    m_x = r_x' M r_x = 125,672.000000 kg
    m_y = r_y' M r_y = 125,672.000000 kg
    m_z = r_z' M r_z = 125,672.000000 kg.

Here each `r` has unit displacement in one translational direction at every node and zero rotations. Summing every translational diagonal indiscriminately would count the physical mass repeatedly; that is not the check used here.

The often-mentioned 0.61 scale would imply 76,659.92 kg for this numerical model. That number is a consequence of the proposed scale, **not an independent mass estimate**. Nor is it justified simply to delete the 2 m overhang contribution: those overhangs are part of the documented structure. The arithmetic `6830*16.4 = 112012 kg` is only an extent sensitivity, not an approved correction.

## 3. Verified train weighing

The original campaign CSVs are included under `data/weighing`. Their wheel masses were checked against Tables 12 and 13 of `sensors_weights_protocols.pdf`, page 11. Each row combines the right and left wheel masses at one axle station; consecutive axle pairs form the eight bogies. The conversion uses the declared convention g = 9.80665 m/s^2.

| Campaign | Sum of weighed mass (kg) | Total static weight (kN) | Mean per bogie (kN) |
|---|---:|---:|---:|
| 2022 | 225,041 | 2,206.898323 | 275.862290 |
| 2025 | 224,071 | 2,197.385872 | 274.673234 |

The exact eight-bogie values are in `verified_static_bogie_loads.csv`. These are campaign **static weighing** data; they are not measured time-varying wheel/rail contact forces during a passage. Using them as prescribed passage forces introduces that explicit modelling assumption.

Across the uploaded results, the median inferred total is 1,647.346274 kN for 2022 and 1,629.317542 kN for 2025. Across all 288 bogie posterior means, the median is 181.259034 kN and 164 values are at or below 182 kN, close to the 180 kN lower bound. This is evidence of tension between the reference-strain inverse and the independent static weights. It does not uniquely identify whether the error is in load representation, observation recovery, structural response, or their combination.

Run 2025/40 has an inferred total of 1,668.866819 kN. Its weighed/inferred total ratio is about 1.3167. Run 2025/22 has a different total ratio, about 1.4006. There is therefore no single demonstrated 1.32 campaign multiplier.

## 4. Why 0.46 -> 0.61 is not yet a computed corrected prediction

The earlier estimate multiplied an observed displacement ratio by a total-load ratio. It is useful as an illustrative scale check. It is not the result of running the current model with the weighed bogie pattern or reconditioning the probabilistic model.

For a linear fixed model, the eight unit-bogie responses combine as

    u(t) = B_u(t) p,       strain(t) = B_epsilon(t) p.

Multiplication of the entire response by a single factor is exact only when every force component receives the same factor and the model and initial conditions remain fixed. The weighed eight-bogie pattern is not generally a common multiple of the inferred pattern. The current conditional prediction also contains a reference-residual correction and load-dependent discrepancy statistics. Those components cannot be corrected by multiplying only the final displacement curve.

The new optional `run_bridgewest_known_load_diagnostic` recomputes the **nominal** unit-bogie basis at the saved passage speed and applies the inferred and weighed bogie vectors separately. It saves both strain and displacement. It does not read PLS/target strains, does not change production predictions, does not re-estimate the posterior, and does not claim calibrated uncertainty. Native execution of this diagnostic has not been performed in the delivery environment.

## 5. Correct interpretation of common structural scaling

For a linear system with unchanged force mapping and zero initial conditions,

    M u_ddot + C u_dot + K u = F,

scaling **all three** matrices by the same positive factor s preserves the eigenvalues and damping ratios and gives `u_s = u/s` under the same F. With a fixed strain observation operator, strains also scale by 1/s. The transformation is not a displacement-only correction.

This mathematical scaling is a possible sensitivity model, not proof that the real bridge has uniformly mis-scaled mass and stiffness. A change to ballast mass alone is not a common scaling of every steel mass term, support stiffness and damping term. Independent total mass also does not uniquely determine its spatial distribution or every stiffness component.

The current shaker branch substitutes the first bending and torsional frequencies/damping in the modal response model while retaining the saved FE shapes/mass. These observed modal quantities do not supply a full physical mass/stiffness identification. Calling the current configuration a detailed, fully physically calibrated bridge model would be too strong.

If loads remain freely inferred from the same reference strains, a common structural scale may be absorbed by a compensating load change. Conversely, imposing static weighing can expose a worse reference fit rather than solve the error. Any accepted physical update must therefore be propagated through **both** strain and displacement, followed by a fresh reference-based inversion and separate target/PLS assessment.

## 6. Evidence that would permit the next physical update

A supported inventory needs the relevant steel section/thickness schedule and member lengths, deck thickness and elevation, net ballast cross-section/extent and density, rail sections and lengths, sleeper type/count/mass, and the walkway/handrail/ancillary contributions. Components already represented by structural elements must not be counted again as added mass. Original bridge-specific CDB/design drawings would be the most direct sources.

An independently force-calibrated shaker frequency-response amplitude could also constrain modal mass/absolute response scale, but that requires an identified force-transfer/input convention and a new amplitude-based calibration. The present frequency/damping state ensemble is not that calibration.

Until such evidence is supplied, the correct delivered value of `bridge_mass_scale_applied` is **1.0**, with `independent_bridge_mass_established = false`. The 0.61 hypothesis is neither accepted nor reported as a physical mass measurement.

## 7. Plotting correction and execution checks

The v3.1 `plot_bridgewest_stabil_results` wrapper delegated only to the displacement plotter. The repaired function defaults to strain and accepts individual events, multiple events, complete campaigns, and campaign/run selections. It uses the actual saved `model_channel_order`, frozen time grids, and assessment-only cached target measurements. A non-ok A070 observation never becomes a primary score merely because a curve is displayed.

Python recomputation against the uploaded native arrays checked 36 passages and 144 target histories, including the 132 quality-qualified histories. Maximum differences from the supplied CSV scores were approximately 5.0e-16 for NRMSE, 2.0e-15 for correlation and 5.0e-16 for coverage. All 36 prediction hashes were unchanged. The details and source hashes are in `verification.json` and `prediction_sha256_before_after.json`.

MATLAB and MATLAB graphics were not executed here. Native MATLAB data/selection regression tests are supplied under `tests`; the Python audit is a check of the inputs, mappings and calculations, not a claim of MATLAB runtime parity. No new 36-passage FE or posterior calculation is claimed.

## Sources

Lorenzen, S. R. (2023). *Railway Bridge Monitoring with Minimal Sensor Deployment: Virtual Sensing and Resonance Curve-Based Drive-by Monitoring*. Chapter 7, Figure 7.2 and Section 7.2.3; Table 5.1 for the explicitly reused steel properties. DOI: 10.26083/tuprints-00026426.

Original dataset metadata: `sensors_weights_protocols.pdf`, page 11, Tables 12–13; `weights_ICE_TD_2022.csv` and `weights_ICE_TD_2025.csv`. Dataset DOI: 10.5281/zenodo.19046962. Original metadata attribution and CC BY-NC-SA 4.0 reuse conditions remain applicable.

User-supplied `sofi-service-4-zekiss-main(1).zip`, inspected for bridge-specific geometry/model inputs; and native `BridgeWest_Stabil_Shaker_v3_results.zip`, used for numerical audit only. Previously generated assumptions/configurations are not independent evidence of physical bridge mass.
