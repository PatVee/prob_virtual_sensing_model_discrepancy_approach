# SOFiSTiK CDB adapter

The supplied `sofi-service-4-zekiss` archive is a generic SOFiSTiK 2020 CDB
interface/service, not the Garather Bach bridge model.  It contains no `.cdb` or
bridge-specific `.dat` file.  It is nevertheless useful because its generated
`sofistik_daten.py` exposes the native records needed to remove assumptions from
the reconstructed Stabil model.

The most relevant records are:

- `20/00` — nodes and support bit codes;
- `100/00` — beam connectivity and property references;
- `200/00` — QUAD connectivity, materials, thicknesses and bedding values;
- `9/NR` — section area, inertias, torsion and material reference;
- `1/NR` — steel material constants and density/nominal weight;
- `30/NR:10` — elastic point springs, direction, stiffnesses and reference area.

`inspect_bridgewest_cdb.py` is an adapter for use on a Windows machine with the
SOFiSTiK 2020 CDB interface DLLs.  Put it beside the supplied `sofistik_connect.py`
and `sofistik_daten.py`, then run, for example:

```text
python inspect_bridgewest_cdb.py Garather_Bach.cdb --out bridgewest_cdb_extract.json
```

No such Garather Bach CDB is present in the current archive, so this adapter has
not been run against the bridge model.  Do not treat its output schema as proof
of a property until the actual CDB has been inspected and its unit system checked.

If the original CDB becomes available, its extracted QUAD thicknesses, beam
sections, support codes and spring records should replace the corresponding
provisional entries in `config/model.json`.  The shaker-informed modal ensemble
then remains unchanged in concept; it is simply regenerated on the better
physical FE model.
