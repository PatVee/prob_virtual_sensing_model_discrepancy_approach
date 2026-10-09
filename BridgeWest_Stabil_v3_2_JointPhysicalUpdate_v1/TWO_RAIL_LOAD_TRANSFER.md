# Bridge West two-wheel / two-rail force-to-response map

This package changes only the **moving-load force application operator**. The FE model, shaker-state ensemble, bounded eight-bogie posterior, A031/B053 conditioning, withheld targets, noise assumptions and posterior bounds remain unchanged.

For bogie resultant \(P_b\), the load chain is

\[
P_b \rightarrow \frac{P_b}{2}\;\text{per axle}
\rightarrow \frac{P_b}{4}\;\text{per wheel}
\rightarrow \text{left/right rail footprints}
\rightarrow \mathbf F_b(t)
\rightarrow \boldsymbol\Phi^T\mathbf F_b(t)
\rightarrow \mathbf q_b(t)
\rightarrow \{\boldsymbol\varepsilon_b(t),\mathbf d_b(t)\}.
\]

The rail ordinates are

\[
y_L=y_c-\frac{g}{2},\qquad y_R=y_c+\frac{g}{2},
\]

with `rail_gauge_m = 1.435`. A uniform rectangular influence kernel is integrated over the tributary rectangles of the deck-shell nodes. Default half dimensions are `2.4 m` longitudinally and `0.6 m` transversely around each rail. These footprint dimensions are **explicit modelling assumptions**, not quantities identified from Bridge West target responses.

With `conserve_wheel_resultant = true`, the footprint weights are renormalised so that every active wheel contributes exactly \(P_b/4\). Thus the operator modifies spatial load transfer without changing the bogie resultant or the original bearing-to-bearing entry timing.

## Files changed

- `shaker_informed/matlab/bw_stabil_force_two_rail.m` — explicit wheel/rail footprint map.
- `shaker_informed/matlab/bw_stabil_force_map.m` — dispatcher retaining a legacy-centreline option.
- `shaker_informed/matlab/bw_stabil_sim_shaker.m` — calls the dispatcher.
- `shaker_informed/matlab/bw_stabil_basis_shaker.m` — records the active load basis.
- `shaker_informed/discrepancy_config.json` — activates and documents the new operator.

To return to the previous operator, set

```json
"load_transfer_model": "centerline_point"
```

No static-weighing information is used by the new force map or posterior.
