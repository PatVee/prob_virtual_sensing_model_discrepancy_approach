> Historical v3/v3.1 document retained for provenance. For the current update, read `V3_2_README.md` and `audit/MASS_AND_LOAD_AUDIT.md`. No mass correction is activated in v3.2.

# Bridge West FE rebuild: findings and limits

## Outcome

A new shell-and-beam Stabil model has been prepared, together with a Python
numerical counterpart and separate plotting scripts. It uses the structural
representation described by Lorenzen for Garather Bach. The available source
material is insufficient to reproduce the original SOFiSTiK model exactly.
The current nominal reconstruction does not reproduce the five modal
observations reliably, so it has not replaced the previous working FE model.

## Correct modal targets

The five observations appear in **Table 7.2, printed page 110**. They were
identified from ambient measurements collected on 3 March 2022, using SSI.
They are historical measured modal properties, not the first five eigenvalues
of a supplied numerical model and not a September 2024/January 2025 shaker fit.

| Label | Measured type | Frequency (Hz) | Damping (%) | MPC |
|---|---|---:|---:|---:|
| 1 | First bending | 6.73 | 1.40 | 0.98 |
| 2 | First torsion | 8.99 | 0.90 | 0.98 |
| 3 | Second bending | 16.14 | 3.60 | 0.84 |
| 4 | Second torsion | 18.68 | 3.00 | 0.93 |
| 5 | Undefined | 24.13 | 3.20 | 0.54 |

MPC is modal phase collinearity; it is not the MAC between an FE shape and a
measured shape. The undefined fifth observation is retained as a diagnostic
target, with its stated quality measure. It is not assigned to a numerical
mode merely because a frequency happens to be nearby. The thesis compares
four identified shapes with FE in Figure 7.14 and treats six selected FE modes
in Figure 7.13. These different counts should not be equated.

Numeric March 2022 ambient mode vectors are not included in the supplied
material. Their MAC values cannot be recomputed from Table 7.2. The printed
MAC values in Figure 7.14 belong to Lorenzen's original model, not this rebuild.

## Material, bearings and track restraint

Section 7.2.3 explicitly refers to the steel properties of Table 5.1: E=210 GPa,
Poisson ratio 0.30, and specific weight 78,500 N/m3. The implementation stores
specific weight as stated and converts it with g=9.80665 m/s2, giving mass
density 8,004.77 kg/m3. This is intentionally distinguished from assuming a
7,850 kg/m3 steel density. The aggregate bridge mass remains a separate nominal
input and includes track/ballast and ancillary components.

The concrete row and the connection springs in Table 5.2 belong to the HUMVIB
footbridge. They are not imported into Bridge West. For Bridge West, the
documented spring quantity is the upper-flange foundation modulus of
700,000 kN/m3 in Section 7.2.3, representing restraint from track/ballast.

The reconstruction retains free rotations at all bearings and fixes their four
vertical translations. A minimal pin/guided/free in-plane restraint arrangement
removes rigid-body motion. Its exact directions remain a support idealisation
until a bearing drawing is available. It is not a claim that the physical
bridge has the identical bearing schedule.

The upper-flange foundation is assembled separately. The nominal case uses
vertical springs over four 1.00 m by 0.50 m flange patches on the overhangs.
Here 0.50 m is an assumed flange width, not a dimension given by the metadata:
the metadata's phrase “0.5 flange width” describes a gauge position.
With the assumed total patch area of 2.00 m2, the summed spring stiffness is
1.40 GN/m. A foundation modulus in N/m3 is multiplied by tributary area; it
is never entered directly as N/m or Nm/rad. The thesis does not specify the
patch dimensions or spring direction sufficiently to claim that this is the
original spring assembly. Alternative patch and direction cases are supplied.

## Geometry and mass

The documented geometry includes the 16.40 m bearing span, the 1.00 m extensions
at either end, and 3.30 m spacing between DT girders from the metadata. Plate
thicknesses, web depth, flange width, deck elevation and cross-member/walkway
sections are provisional nominal parameters. The full parameter/source mapping
is in `evidence/parameter_provenance.csv`.

The coarse model has 963 nodes, 864 shell elements, 188 beam elements and 5,771
free DOFs. It is not intended to reproduce the original 585,084-DOF mesh by
matching its element count. The nominal total mass is 125,672 kg over the full
18.40 m model: approximately 24,427 kg steel, 883 kg walkway surface mass, and
100,361 kg residual track/ballast mass. These values depend on the provisional
sections and the assumption that 6,830 kg/m applies over the modelled length.

## Modal comparison and refinement

| Thesis observation | Nominal FE candidate | Difference | Assessment |
|---|---:|---:|---|
| First bending, 6.73 Hz | 5.21 Hz | −22.57% | Frequency does not match |
| First torsion, 8.99 Hz | 7.74 Hz | −13.94% | Frequency does not match |
| Second bending, 16.14 Hz | 23.72 Hz | +46.99% | Candidate mapping is unstable |
| Second torsion, 18.68 Hz | Unassigned | — | No reliable candidate |
| Undefined, 24.13 Hz | Unassigned | — | No forced correspondence |

Mode screening uses the relative motion of the two girders, vertical response,
upper-flange amplitude and longitudinal lobe count. It is a diagnostic screen,
not a replacement for experimental MAC. The nominal candidate indices are 3,
8 and 109; local/lateral modes account for the intervening eigenvalues.

On halving the mesh-size limits, the first bending candidate changes from
5.2110 to 5.1461 Hz (−1.25%) and the first torsion candidate from 7.7364 to
7.6849 Hz (−0.67%). Their common-node translational MAC values are 0.9996
and 0.9643. These two shapes are stable in this refinement.

The coarse second-bending candidate has a best common-node match at 32.2872 Hz
on the finer mesh, with MAC only 0.7561. Its higher-mode identification is
therefore unresolved. A similar upper-girder profile alone can select a different
local/deck mode; full common-node motion was used for the refinement check.
These results do not justify fitting all five frequencies by freely adjusting
unknown geometry or by changing eigenvalues independently of K and M.

An additional comparison against the supplied 2025 four-station shaker shapes
gives MAC 0.9987 and 0.9961 for the first bending and torsional candidates.
This is a different measurement campaign from the thesis. High MAC for those
few coordinates does not establish correct frequency, absolute stiffness,
response amplitude or higher-mode behaviour.

Removing the upper-flange restraint changes the first bending/torsional
candidates to approximately 4.92/7.27 Hz. The direction and patch sensitivity
results are in `assumption_and_mesh_sensitivity.csv`; candidate labels in those
cases remain screening outputs. A frequency-only match in a sensitivity case
must not be accepted without shape and geometry checks.

## Numerical verification

The Python element implementation passed constant-curvature shell energy and
analytical simply supported plate checks. The plate first-frequency error
reduced from 0.75% on the 4-by-4 mesh to 0.11% on the 12-by-12 mesh. A separate
cantilever-beam frequency agreed with the analytical value to better than
0.0001%. These checks validate basic mechanics, not the physical bridge model.

The assembled model conserves total translational mass, has symmetric K/M,
zero prescribed bearing displacement, vertical static equilibrium error below
2e-9 N under the 1 kN diagnostic load, and maximum relative eigenpair residual
about 2.6e-9. No artificial drilling mass was introduced. See the JSON checks
for unrounded values.

MATLAB syntax was checked. MATLAB/Octave execution was unavailable; the numerical
results come from Python. `check_bridgewest_lorenzen_parity` is supplied for
checking the assembled MATLAB matrices on the user's machine. MATLAB runtime
parity is not yet an established result.

## Completion requirements

The next useful input is the original SOFiSTiK `.dat`/`.cdb` model or dimensioned
bridge drawings. The required details are the web/flange/deck thickness schedule,
web depth and deck elevation, cross-beam and stiffener dimensions and connections,
walkway/handrail sections, bearing restraint directions, ballast distribution,
and upper-flange spring direction and footprint. Numerical ambient mode vectors
would enable the four literature mode assignments to be verified by MAC.

After these are supplied, update the geometry, establish mesh-stable modal
correspondence, and only then regenerate the forward simulations for the
existing fixed-library probabilistic method. The current package does not claim
improved B081 or displacement prediction and does not alter that methodology.
