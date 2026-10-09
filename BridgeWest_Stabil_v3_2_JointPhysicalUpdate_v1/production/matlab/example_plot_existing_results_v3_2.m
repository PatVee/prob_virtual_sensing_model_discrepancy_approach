% Edit paths; this script plots saved results and does NOT rerun inference.
outputRoot = 'D:\BridgeWest_Stabil_Shaker_v3_results';
root=fileparts(mfilename('fullpath'));addpath(root,'-begin');
plot_bridgewest_stabil_results(outputRoot,'Event','2025_run22');
plot_bridgewest_stabil_results(outputRoot,'Event','2022_run20');
% All saved passages:
% plot_bridgewest_stabil_results(outputRoot,'Event','all','Visible','off','CloseAfterSave',true);
