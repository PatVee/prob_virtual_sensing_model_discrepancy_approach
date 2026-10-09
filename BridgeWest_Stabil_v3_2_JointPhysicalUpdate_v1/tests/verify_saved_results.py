"""Recompute the v3.2 audit from native saved MATLAB arrays; never edits results.

Usage:
    python tests/verify_saved_results.py --results-root PATH_TO_RESULTS
Requires numpy, pandas, scipy, h5py. This does not execute MATLAB graphics.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.sparse import csc_matrix

LABELS = ('A058', 'A070', 'B031', 'B081')
CHANNELS = ('str_x_BA_uf_058', 'str_x_BA_uf_070', 'str_x_BB_uf_031', 'str_x_BB_uf_081')
OFFSETS = np.array([0, 2.6, 19, 21.6, 25.9, 28.5, 44.9, 47.5, 51.8, 54.4, 70.8, 73.4, 77.7, 80.3, 96.7, 99.3])
G = 9.80665


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def strings(f: h5py.File, key: str) -> list[str]:
    return [''.join(chr(int(c)) for c in f[ref][()].ravel())
            for ref in f[key][()].ravel()]


def verify(results: Path, package: Path, audit: Path) -> dict:
    audit.mkdir(parents=True, exist_ok=True)
    files = sorted((results / 'predictions').glob('*_run*.mat'))
    assert len(files) == 36, f'Expected the supplied 36 passages, found {len(files)}'
    events = [f.stem for f in files]
    assert sum(x.startswith('2022') for x in events) == 12
    assert sum(x.startswith('2025') for x in events) == 24
    assert '2025_run22' in events and '2022_run20' in events
    assert '2022_run22' not in events
    before = {f.name: sha256(f) for f in files}
    raw_scores = pd.read_csv(results / 'tables' / 'target_assessment_all.csv')
    weights, weigh_rows, bogie_rows = {}, [], []
    sources = []
    for year in (2022, 2025):
        source = package / 'data' / 'weighing' / f'weights_ICE_TD_{year}.csv'
        w = pd.read_csv(source).sort_values('no')
        assert np.array_equal(w['no'].values, np.arange(1, 17))
        assert np.allclose(w.location_m, OFFSETS, atol=1e-10, rtol=0)
        axle = w.weight_right_kg.values + w.weight_left_kg.values
        p = axle.reshape(8, 2).sum(axis=1) * G / 1000
        weights[year] = p
        weigh_rows.append(dict(campaign=year, total_mass_kg=int(axle.sum()),
                               gravity_m_s2=G, total_load_kN=float(p.sum()),
                               mean_bogie_load_kN=float(p.mean())))
        for i, pi in enumerate(p, 1):
            bogie_rows.append(dict(campaign=year, bogie=i, static_load_kN=float(pi)))
        sources.append(dict(source=str(source.relative_to(package)), sha256=sha256(source),
                            role='original independent static campaign weighing'))
    metric_rows, load_rows, loads_all = [], [], []
    max_errors = np.zeros(3)
    for path in files:
        event = path.stem
        with h5py.File(path, 'r') as f:
            t = f['time_s'][()].ravel()
            assert np.isfinite(t).all() and (np.diff(t) > 0).all()
            order = strings(f, 'model_channel_order')
            idx = [order.index(ch) for ch in LABELS]
            mean = f['strain_mean_microstrain'][()].T[:, idx]
            lo = f['strain_lo_microstrain'][()].T[:, idx]
            hi = f['strain_hi_microstrain'][()].T[:, idx]
            assert mean.shape == lo.shape == hi.shape == (len(t), 4)
            assert all(np.isfinite(x).all() for x in (mean, lo, hi))
            assert (lo <= hi).all()
            assert not bool(f['target_measurements_loaded_during_inference'][()].item())
            p = f['load_mean_kN'][()].ravel()
            loads_all.extend(p)
            year = int(event[:4]); total = float(p.sum()); wp = float(weights[year].sum())
            load_rows.append(dict(event=event, campaign=year,
                speed_kmh=float(f['speed_m_s'][()].item()) * 3.6,
                inferred_total_kN=total, weighed_total_kN=wp,
                inferred_over_weighed=total / wp, weighed_over_inferred=wp / total,
                n_bogies_le182_kN=int((p <= 182).sum())))
        z = loadmat(results / 'validation_targets' / f'{event}.mat', simplify_cells=True)
        assert len(t) == len(z['time']) and np.max(np.abs(t - z['time'])) < 1e-8
        assert tuple(z['channels']) == CHANNELS
        for j, label in enumerate(LABELS):
            y = z['Y'][:, j]; flag = str(z['flags'][j]); qualified = bool(flag == 'ok' and np.isfinite(y).all())
            nrmse = float(np.sqrt(np.mean((mean[:, j] - y)**2)) / max(np.ptp(y), np.finfo(float).eps))
            correlation = float(np.corrcoef(mean[:, j], y)[0, 1])
            coverage = float(np.mean((y >= lo[:, j]) & (y <= hi[:, j])))
            q = raw_scores.loc[(raw_scores.event == event) & (raw_scores.target == label)]
            assert len(q) == 1 and qualified == bool(q.iloc[0].primary_score)
            ref = q[['nrmse', 'correlation', 'coverage']].iloc[0].values
            err = np.abs(np.array([nrmse, correlation, coverage]) - ref)
            assert np.max(err) < 1e-10, (event, label, err)
            max_errors = np.maximum(max_errors, err)
            metric_rows.append(dict(event=event, target=label, quality_flag=flag,
                primary_score=qualified, nrmse=nrmse, correlation=correlation, coverage=coverage))
    fe_path = results / 'fe_model' / 'bridgewest_lorenzen_fe.mat'
    with h5py.File(fe_path, 'r') as f:
        fe = f['fe']; m = fe['M']; n = int(m.attrs['MATLAB_sparse'])
        M = csc_matrix((m['data'][()].ravel(), m['ir'][()].ravel(), m['jc'][()].ravel()), shape=(n, n))
        masses = []
        for axis in range(3):
            r = np.zeros(n); r[axis::6] = 1
            masses.append(float(r @ M @ r))
        component_values = [float(fe[x][()].item()) for x in ('steel_mass', 'walkway_mass', 'residual_mass')]
        nominal = float(fe['c']['mass_per_length_total_kg_m'][()].item()) * float(fe['c']['length_m'][()].item())
    assert np.allclose(masses, nominal, atol=1e-6, rtol=0)
    assert abs(sum(component_values) - nominal) < 1e-6
    after = {f.name: sha256(f) for f in files}
    assert before == after, 'Protected prediction files changed!'
    metrics = pd.DataFrame(metric_rows); loads = pd.DataFrame(load_rows)
    assert len(metrics) == 144 and metrics.primary_score.sum() == 132
    metrics.to_csv(audit / 'recomputed_strain_scores.csv', index=False)
    loads.to_csv(audit / 'train_load_audit.csv', index=False)
    pd.DataFrame(weigh_rows).to_csv(audit / 'verified_train_weighing_totals.csv', index=False)
    pd.DataFrame(bogie_rows).to_csv(audit / 'verified_static_bogie_loads.csv', index=False)
    mass_table = pd.DataFrame(dict(component=['steel', 'walkway_surface', 'residual_added_mass'],
                                  mass_kg=component_values,
                                  fraction_of_total=np.array(component_values) / nominal))
    mass_table.to_csv(audit / 'bridge_mass_components.csv', index=False)
    sources.append(dict(source='uploaded results/fe_model/bridgewest_lorenzen_fe.mat',
                        sha256=sha256(fe_path), role='saved numerical model, not independent physical mass'))
    result = dict(
        executed_in='Python 3 using native MATLAB v7/v7.3 arrays',
        matlab_execution_performed=False, matlab_graphics_execution_performed=False,
        full_FE_campaign_rerun_performed=False,
        checked_passages=len(events), passages_2022=12, passages_2025=24,
        checked_target_histories=len(metrics), primary_quality_qualified_histories=int(metrics.primary_score.sum()),
        maximum_score_recomputation_errors=dict(zip(['nrmse','correlation','coverage'], max_errors.tolist())),
        rigid_translation_mass_kg=dict(zip(['x','y','z'], masses)),
        imposed_total_mass_kg=nominal,
        mass_partition_sum_kg=sum(component_values),
        independent_bridge_mass_established=False,
        bridge_mass_scale_applied=1.0,
        all_288_bogie_posterior_means_median_kN=float(np.median(loads_all)),
        bogie_means_le182_kN=int(np.sum(np.array(loads_all) <= 182)),
        campaign_median_load_totals_kN={str(year): float(loads.loc[loads.campaign == year, 'inferred_total_kN'].median()) for year in (2022, 2025)},
        prediction_files_unchanged=before == after,
        source_records=sources)
    (audit / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
    (audit / 'prediction_sha256_before_after.json').write_text(json.dumps(before, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-root', required=True, type=Path)
    parser.add_argument('--package-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--audit-dir', type=Path)
    args = parser.parse_args()
    out = args.audit_dir or args.package_root / 'audit'
    print(json.dumps(verify(args.results_root, args.package_root, out), indent=2))
