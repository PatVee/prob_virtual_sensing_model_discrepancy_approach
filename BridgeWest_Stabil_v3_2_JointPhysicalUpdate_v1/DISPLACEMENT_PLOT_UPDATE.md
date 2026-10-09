# Displacement plot update

The production workflow now carries the A031/B053-conditioned posterior to vertical displacement as well as withheld strain.

For every successfully predicted passage, the prediction MAT file stores:

- BA mean, standard deviation and nominal 95% interval at 41 longitudinal stations;
- BB mean, standard deviation and nominal 95% interval at the same 41 stations;
- the exact dataset-coordinate vector obtained from the detailed Stabil observation map.

The displacement measurement remains completely outside inference. For 2025 passages with a PLS record, `score_bridgewest_stabil_displacement.m` opens the PLS file only after the prediction exists and creates an independent measured own-midspan-peak spatial profile for BA and BB. It then reports spatial RMSE, range-normalised RMSE, spatial correlation and the predicted/measured peak-magnitude ratio.

`plot_bridgewest_stabil_displacement.m` produces two outputs:

1. a selected-event BA/BB spatial displacement figure showing the predicted own-peak profile, nominal 95% interval and measured PLS own-peak profile when available;
2. a 2025 campaign summary showing peak-amplitude transfer and spatial NRMSE for BA and BB.

The default `run_all_bridgewest_stabil_v3.m` workflow now executes strain assessment, displacement assessment and displacement plotting after the reference-only prediction stage.
