from pathlib import Path
base=Path(__file__).resolve().parent
m=base/'production'/'matlab'
checks={
 'plot_function_exists': (m/'plot_bridgewest_stabil_displacement.m').is_file(),
 'plot_wrapper_exists': (m/'plot_bridgewest_stabil_results.m').is_file(),
 'displacement_scorer_exists': (m/'score_bridgewest_stabil_displacement.m').is_file(),
 'pls_reader_exists': (m/'raw_functions'/'bw_pls_raw_snapshot.m').is_file(),
}
runner=(m/'run_bridgewest_stabil_shaker_informed.m').read_text()
runall=(m/'run_all_bridgewest_stabil_v3.m').read_text()
scorer=(m/'score_bridgewest_stabil_displacement.m').read_text()
checks.update({
 'prediction_saves_coordinates': 'displacement_x_dataset_m' in runner,
 'prediction_keeps_BA_BB_fields': all(x in runner for x in ['disp_BA_mean_mm','disp_BB_mean_mm']),
 'runall_calls_displacement_assessment': 'score_bridgewest_stabil_displacement' in runall,
 'runall_calls_plotting': 'plot_bridgewest_stabil_results' in runall,
 'pls_is_post_prediction': 'PLS_first_opened_after_prediction' in scorer,
 'pls_not_used_for_inference': 'prediction_stage_uses_PLS' in scorer,
})
for k,v in checks.items(): print(f'{k}: {v}')
if not all(checks.values()): raise SystemExit('Displacement update audit failed.')
print('All displacement update static checks passed.')
