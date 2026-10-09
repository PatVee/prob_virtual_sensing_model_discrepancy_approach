"""Bridge West SOFiSTiK CDB extraction adapter.

This adapter is intended to be placed beside the user-supplied `sofi-service`
Python files (`sofistik_connect.py`, `sofistik_daten.py`) on a Windows machine
with the SOFiSTiK 2020 CDB interface DLLs available.  It is not executed in the
current package because no Garather Bach `.cdb` file is supplied.

The objective is provenance recovery, not recalibration: export native nodes,
QUAD thicknesses/materials, beam property references and elastic point springs
so that provisional Stabil assumptions can be replaced by source-model values.
"""
from __future__ import annotations
from ctypes import byref, c_int, sizeof
from pathlib import Path
import argparse, json

from sofistik_connect import connect_to_cdb, close_cdb, py_sof_cdb_get
from sofistik_daten import (
    CNODE, cnode, CQUAD, cquad, CBEAM, cbeam, CGPT_SPR, cgpt_spr,
    CSECT, csect, CMAT_STEE, cmat_stee,
)


def read_records(index, kwh, kwl, obj, cls, convert):
    rows=[]; ie=c_int(0); reclen=c_int(sizeof(cls))
    while ie.value < 2:
        reclen.value=sizeof(cls)
        ie.value=py_sof_cdb_get(index, int(kwh), int(kwl), byref(obj), byref(reclen), 1)
        if ie.value < 2:
            row=convert(obj)
            if row is not None: rows.append(row)
    return rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('cdb', type=Path)
    ap.add_argument('--out', type=Path, default=Path('bridgewest_cdb_extract.json'))
    args=ap.parse_args()
    if not args.cdb.exists(): raise FileNotFoundError(args.cdb)
    idx,stat=connect_to_cdb(str(args.cdb))
    if stat.value==0: raise RuntimeError('SOFiSTiK CDB did not open')
    try:
        nodes=read_records(idx,20,0,cnode,CNODE,lambda r:{
            'nr':int(r.m_nr),'kfix':int(r.m_kfix),'xyz':[float(x) for x in r.m_xyz]})
        quads=read_records(idx,200,0,cquad,CQUAD,lambda r:{
            'nr':int(r.m_nr),'nodes':[int(x) for x in r.m_node],'material':int(r.m_mat),
            'type':int(r.m_nra),'thickness':[float(x) for x in r.m_thick],
            'bedding_normal':float(r.m_cb),'bedding_tangent':float(r.m_cq)})
        beams=read_records(idx,100,0,cbeam,CBEAM,lambda r:{
            'nr':int(r.m_nr),'nodes':[int(x) for x in r.m_node],
            'property_number':int(r.m_np),'length':float(r.m_dl)})
        # The supplied SOFiSTiK 2020 schema labels 30/NR:10 as elastic spring support.
        springs=read_records(idx,30,10,cgpt_spr,CGPT_SPR,lambda r:{
            'nr':int(r.m_nr),'group':int(r.m_ng),'node':int(r.m_ipc),'direction_code':int(r.m_iri),
            'direction':[float(x) for x in r.m_t],'reference_area':float(r.m_aref),
            'cp':float(r.m_cp),'cq':float(r.m_cq),'cm':float(r.m_cm),'gap':float(r.m_gap)})
        prop_numbers=sorted({b['property_number'] for b in beams if b['property_number']>0})
        sections=[]
        for nr in prop_numbers:
            rec=read_records(idx,9,nr,csect,CSECT,lambda r:{
                'nr':nr,'id':int(r.m_id),'material':int(r.m_mno),'area':float(r.m_a),
                'Ay':float(r.m_ay),'Az':float(r.m_az),'It':float(r.m_it),
                'Iy':float(r.m_iy),'Iz':float(r.m_iz),'Iyz':float(r.m_iyz),
                'E':float(r.m_em),'G':float(r.m_gm),'nominal_weight':float(r.m_gam)})
            sections.extend(rec)
        mat_numbers=sorted({q['material'] for q in quads if q['material']>0} | {s['material'] for s in sections if s['material']>0})
        materials=[]
        for nr in mat_numbers:
            rec=read_records(idx,1,nr,cmat_stee,CMAT_STEE,lambda r:{
                'nr':nr,'id':int(r.m_id),'type':int(r.m_type),'E':float(r.m_emod),
                'nu':float(r.m_mue),'G':float(r.m_gmod),'nominal_weight':float(r.m_gam),
                'density':float(r.m_rho)})
            materials.extend(rec)
        out={'source_cdb':str(args.cdb.resolve()),'nodes':nodes,'quads':quads,'beams':beams,
             'elastic_point_springs':springs,'sections':sections,'steel_materials':materials,
             'note':'Verify SOFiSTiK unit settings and CDB key variants against the native model before importing into Stabil.'}
        args.out.write_text(json.dumps(out,indent=2))
        print(f"wrote {args.out}: {len(nodes)} nodes, {len(quads)} quads, {len(beams)} beams, {len(springs)} springs")
    finally:
        close_cdb(idx,stat)

if __name__=='__main__': main()
