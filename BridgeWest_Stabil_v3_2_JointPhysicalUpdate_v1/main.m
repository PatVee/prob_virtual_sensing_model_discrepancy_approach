% Bridge West literature-based FE rebuild. MATLAB R2021a or newer.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab'));
fe=run_bridgewest_lorenzen();
figs=plot_bridgewest_lorenzen(fe);
