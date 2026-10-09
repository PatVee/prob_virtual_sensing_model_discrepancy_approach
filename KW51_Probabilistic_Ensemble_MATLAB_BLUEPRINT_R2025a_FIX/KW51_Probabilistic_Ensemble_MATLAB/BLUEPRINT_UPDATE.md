# Blueprint-BC update

This revision makes the official KW51 blueprint support arrangement the only
production FE boundary condition used by the probabilistic ensemble.

Key changes:

- bundled `fe_models/model_15.mat` replaced by the nominal blueprint bank
  (`fixed=0`, `kxx=0`, two Leuven longitudinal support coordinates);
- new `build_blueprint_bank.m` replaces the four-bearing floating builder;
- each run copies/rebuilds its FE bank into `run_*/fe_bank/model_15.mat`;
- `kwpe_simulate.m` derives `nSupport` from `B.G`, requires exactly two
  support columns, and does not subtract the support-force mean;
- the separate common-translation prior is disabled (`translation_sd_mm=0`);
- manuscript posterior-load replay uses the run-local bank rather than the
  package-level bank;
- selected-event and standard result plots show the blueprint predictive
  interval directly; the obsolete translation-sensitivity figure is removed;
- axes toolbars are hidden before MATLAB figure export;
- `main.m` no longer loads a stale `out.mat` and the numeric event-selection
  example is corrected;
- Python plotting counterparts follow the same run-local blueprint bank;
- the FE numerical fixture is regenerated for two Leuven support histories.

The package deliberately rejects legacy four-support floating banks/draws
instead of silently mixing them with blueprint results.


## R2025a figure-export compatibility

The manuscript, selected-event and run-time figure exporters now use `exportgraphics` rather than legacy `print(fig,...)`. This avoids the MATLAB R2025a `checkArgsForHandleToPrint / Invalid graphics object` failure observed with tiled layouts. Each exporter verifies the figure handle, hides axes toolbars, writes a 220-dpi PNG, writes a vector PDF where supported, and falls back to an image PDF if vector export is unavailable. This is a plotting-only change and does not alter inference or validation results.
