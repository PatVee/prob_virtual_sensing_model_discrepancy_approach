# Stabil routines

The MATLAB files in `Stabil` were copied, without changes, from the Stabil
dependency supplied in the user's existing KW51 package. Their original authors'
headers and references are retained. No new licence is asserted for those files.
Their SHA-256 hashes are in `evidence/stabil_source_hashes.json`.

The shell4 implementation uses an incompatible-mode bilinear membrane and four
overlaid DKT triangles. The Python implementation is a numerical translation for
the planar elements in this rebuild. It is not a general substitute for Stabil.
The Python beam corresponds to the MATLAB option disabling bending rotary
inertia, with infinite shear factors and retained torsional mass inertia.

The original complete Stabil installation may be substituted for this subset.
