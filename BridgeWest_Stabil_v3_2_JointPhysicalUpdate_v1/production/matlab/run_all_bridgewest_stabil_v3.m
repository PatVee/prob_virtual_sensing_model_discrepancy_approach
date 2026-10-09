
clc; clear all; close all;

root=fileparts(mfilename('fullpath'));addpath(root);
% Edit only these two paths, then run this script in MATLAB.
% dataRoot = 'D:\BridgeWestData';
% outputRoot = 'D:\BridgeWest_Stabil_Shaker_v3_1_results';
dataRoot = 'C:\Users\cvvi\Downloads\BridgeWest_Stabil_v3_2_MassAudit_Plots\BridgeWest_Stabil_v3_2_MassAudit_Plots\data_main';
outputRoot = 'C:\Users\cvvi\Downloads\BridgeWest_Stabil_v3_2_MassAudit_Plots\BridgeWest_Stabil_v3_2_MassAudit_Plots\BridgeWest_Stabil_Shaker_v3_results';

% Stage 1: 36 reference-only predictions. Only A031/B053 are opened here.
run_bridgewest_stabil_shaker_informed(dataRoot,outputRoot,'Campaign','all','ContinueOnError',false);
audit_bridgewest_reference_only_protocol(outputRoot);

% Stage 2: separate post-prediction assessment. Targets and PLS are first opened here.
score_bridgewest_stabil_results(dataRoot,outputRoot);
score_bridgewest_stabil_displacement(dataRoot,outputRoot);
audit_bridgewest_reference_only_protocol(outputRoot);

% Stage 3: plots read saved predictions + saved assessment only.
plot_bridgewest_stabil_results(outputRoot,'DataRoot',dataRoot,'Event','all', ...
    'Visible','on','CloseAfterSave',true);
plot_bridgewest_stabil_results(outputRoot,'Event','2022_run20', ...
    'Figures','displacement');

plot_bridgewest_stabil_results(outputRoot, 'Event','2025_run22', ...
    'Figures',{'strain','displacement'});
plot_bridgewest_stabil_results(outputRoot, ...
    'Event',{'2022_run20','2025_run22'});
plot_bridgewest_stabil_results(outputRoot, ...
    'Campaign',2025, ...
    'Runs',[20 22 40]);