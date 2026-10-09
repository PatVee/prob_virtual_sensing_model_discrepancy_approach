from __future__ import annotations
import json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROD = ROOT / 'production'
MAT = PROD / 'matlab'
CFG = PROD / 'config'

def must(cond, msg):
    if not cond:
        raise AssertionError(msg)

ev = json.loads((CFG/'event_sets.json').read_text())
pc = json.loads((CFG/'prediction_config.json').read_text())
must(len(ev['runs_2022']) == 12, '2022 event count must be 12')
must(len(ev['runs_2025']) == 24, '2025 event count must be 24')
must(len(ev['runs_2022']) + len(ev['runs_2025']) == 36, 'total event count must be 36')
must(pc['measured_reference_channels'] == ['str_x_BA_uf_031','str_x_BB_uf_053'], 'reference channels changed')
must(pc['target_labels'] == ['A058','A070','B031','B081'], 'target labels changed')
must(pc['target_measurements_used_in_inference'] is False, 'target inference flag must be false')
must(pc['additional_target_derived_discrepancy'] is False, 'target-derived discrepancy must be false')

runner = (MAT/'run_bridgewest_stabil_shaker_informed.m').read_text()
prep = (MAT/'raw_functions'/'bw_prepare_raw_event.m').read_text()
score = (MAT/'score_bridgewest_stabil_results.m').read_text()
select = (MAT/'raw_functions'/'bw_readstd_selected.m').read_text()
kin = (MAT/'raw_functions'/'bw_kinematics.m').read_text()

must('bw_prepare_raw_event(A,eid,pre,outputRoot,refsMeasured)' in runner, 'runner must pass refs explicitly')
must("'target_measurements_loaded_during_inference',false" in runner, 'prediction target-use flag missing')
must('bw_prepare_validation_targets' not in runner, 'runner must not open validation targets')
must('bw_prepare_validation_targets' in score, 'score must be the target-opening stage')
must("explicitly permitted inference strains" in prep, 'strict prep header missing')
must('bw_readstd_selected(file,chs)' in (MAT/'raw_functions'/'bw_prepstrain.m').read_text(), 'strain prep must use selective reader')
must('bw_readstd_selected(file,railNames)' in kin, 'kinematics must use selective rail reader')
must("h5read(file,'/data'" in select and "requested" in select, 'selective HDF5 reader missing')

# There must be no production configuration instruction to estimate a second/local discrepancy from targets.
for p in [ROOT/'README.md', ROOT/'PROTOCOL.md', PROD/'README.md', CFG/'prediction_config.json']:
    txt = p.read_text().lower()
    must('estimate the second' not in txt, f'obsolete second-discrepancy wording in {p}')

# Basic MATLAB source sanity: function files should end with an end statement.
for p in list(MAT.glob('*.m')) + list((MAT/'raw_functions').glob('*.m')):
    txt = p.read_text().strip()
    if txt.startswith('function '):
        must(re.search(r'\bend\s*$', txt), f'{p.name} does not end with end')

print(json.dumps({
    'package': ROOT.name,
    'events_2022': 12,
    'events_2025': 24,
    'events_total': 36,
    'reference_channels': ['A031','B053'],
    'targets': ['A058','A070','B031','B081'],
    'strict_selected_hdf5_reader': True,
    'target_opening_stage': 'score_bridgewest_stabil_results only',
    'additional_target_derived_discrepancy': False,
    'static_checks': 'PASS'
}, indent=2))
