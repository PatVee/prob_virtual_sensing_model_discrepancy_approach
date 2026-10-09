# KW51 probabilistic ensemble — blueprint boundary-condition package

This package performs the complete native-MATLAB KW51 workflow from the
original monitoring ZIP archives:

1. prepare permitted reference/rail inputs;
2. infer passage kinematics from the permitted measurements;
3. generate 900 load–stiffness candidates;
4. generate 48 paired quasi-static FE discrepancy simulations;
5. condition the ensemble with the three reference strains;
6. predict the two withheld strains and the two Leuven bearing displacements;
7. validate only after the predictions have been frozen; and
8. produce manuscript and selected-event figures.

## Boundary condition used

The package now uses the **official KW51 blueprint boundary condition** as the
nominal FE model (`fixed=0`, `kxx=0`). The operative support arrangement is:

- Brussels N1: x, y and z restrained;
- Brussels N2: x and z restrained;
- Leuven N3: y and z restrained, longitudinal x free;
- Leuven N4: z restrained, longitudinal x and transverse y free.

The protected displacement channels `dABLA` and `dABLC` are the longitudinal
Leuven bearing coordinates N3-x and N4-x. The model-discrepancy support term
therefore uses **two** equivalent longitudinal support-force histories, one at
each Leuven bearing. Their mean is **not** removed, because the Brussels end
already provides the longitudinal reference.

The previous four-bearing floating model, its zero-mean longitudinal gauge,
and the independent common-translation prior are not used in this package.
For compatibility with old result-table schemas, `translation_sd_mm` is retained
in settings/validation files but is fixed to zero and has no effect.

## Run

Put the original raw archives in `data/` and run:

```matlab
main
```

or call directly:

```matlab
out = run_kw51_probabilistic_matlab(dataSource,outputRoot, ...
    'Events','all', ...
    'Validate',true, ...
    'Plot',true);
```

The package includes a prebuilt blueprint response bank. To regenerate it from
the included official FE/Stabil sources:

```matlab
out = run_kw51_probabilistic_matlab(dataSource,outputRoot, ...
    'Events','all', ...
    'Validate',true, ...
    'Plot',true, ...
    'RebuildFE',true);
```

`TranslationSDmm` is retained only for call compatibility and must be zero.

## FE-bank integrity

Every new run copies or rebuilds the FE bank into:

```text
<run folder>/fe_bank/model_15.mat
```

Inference, prediction replay, posterior-load plotting and manuscript plotting
all use this exact run-local file. The code refuses a bank unless `G` has two
support columns and the bank is explicitly identified as the blueprint BC.
This prevents a two-support blueprint run from being mixed with the previous
four-support floating bank.

## Reference and withheld channels

Reference strains used for conditioning:

- `sgBD1011A`
- `sgBD1415A`
- `sgBD2728A`

Withheld validation outputs:

- `sgBD1718C`
- `sgBD2324C`
- `dABLA`
- `dABLC`

Same-station exclusions remain `sgBD1718A`, `sgBD2324AB` and `sgBD2324AT`.
The protected measurements are not used in passage fitting, candidate
weighting, discrepancy construction or posterior conditioning.

## Ensemble settings

The default calculation retains the established probabilistic design:

- 900 load–stiffness candidates;
- stiffness response scales 0.8, 0.9 and 1.0;
- 48 paired nominal/perturbed FE histories;
- load, kinematic, response-scale, strain-recovery and Leuven support-force
  perturbations;
- reference-only Gaussian likelihood weighting and cross-output conditioning;
- 256 coherent conditional response paths;
- nominal pointwise 95% predictive intervals.

The 800 kN support-force scale is a prescribed discrepancy scale. It is not a
measured bearing force or an identified friction law.

## Plotting

All manuscript figures:

```matlab
plot_kw51_manuscript_figures(out,'Visible','on');
```

One accepted passage, using chronological event number 1, 2, 3, ...:

```matlab
plot_kw51_selected_event(out,1,'Visible','on');
```

A full event ID is also accepted.

The manuscript plotting routine reads the FE bank from the run folder rather
than from the package-level `fe_models` folder. The obsolete common-translation
sensitivity plot is not produced for blueprint-BC runs.

Python counterparts are included:

```python
from plot_kw51_manuscript_figures import plot_kw51_manuscript_figures
from plot_kw51_selected_event import plot_kw51_selected_event

plot_kw51_manuscript_figures("KW51ProbabilisticResults")
plot_kw51_selected_event("KW51ProbabilisticResults", 1)
```

## Outputs

Each execution creates a timestamped run folder containing:

- `fe_bank/`: exact blueprint FE bank used by that run;
- `inputs/`: freshly prepared permitted inputs;
- `predictions/`: posterior prediction products;
- `draws/`: load/stiffness and paired-discrepancy draws;
- `validation/`: withheld-output comparisons;
- `validation_scores.csv`, `validation_summary.csv`, `validation_skipped.csv`;
- `reference_diagnostics.csv`;
- `settings.mat`, `settings.json`;
- prediction-freeze hashes and execution log;
- `figures/` and optional `manuscript_figures/`.

## Requirements and verification

MATLAB R2018b or newer and Signal Processing Toolbox are required. Stabil and
the official KW51 FE sources are included.

This revised package was assembled and statically audited outside MATLAB. A
MATLAB runtime was not available in the preparation environment, so the full
raw-ZIP workflow has not been re-executed here. The bundled blueprint bank is
the nominal `[1,1,1]` response bank produced previously from the same official
`fixed=0`, `kxx=0` formulation; the self-test FE fixture has been regenerated
for its two Leuven support coordinates.

Run the supplied numerical tests before the full study:

```matlab
tests/test_kw51_probabilistic
```

Optionally exercise a fresh kinematic fit as well:

```matlab
test_kw51_probabilistic('RunFit',true)
```

## Attribution

KW51 monitoring data: Maes and Lombaert, DOI 10.5281/zenodo.3745914.
KW51 FE model: Maes, van de Velde and Lombaert, DOI 10.5281/zenodo.17366772.
Stabil: François et al., DOI 10.1002/cae.22391.
Published virtual-sensing validation: Maes and Lombaert,
DOI 10.1016/j.ymssp.2023.110142.
