% Run this script, replacing your earlier main/comparison call.
clc;clear;close all;
zipRoot='C:\Users\cvvi\Downloads\Objective2-VituralS-disp\KW51_Probabilistic_Ensemble_MATLAB_BLUEPRINT_R2025a_FIX';
addOnRoot=fileparts(mfilename('fullpath'));addpath(addOnRoot,'-begin');
originalRoot=fullfile(zipRoot,'KW51_Probabilistic_Ensemble_MATLAB');
dataRoot=fullfile(originalRoot,'data');
passageZIPs={fullfile(dataRoot,'traindata_201912_part1.zip'); ...
             fullfile(dataRoot,'traindata_201912_part2.zip')};
% If your files end in _copy.zip, change only their names above.
outputRoot=fullfile(zipRoot,'KW51_modal_only_results');
out=run_kw51_modal_only(originalRoot,passageZIPs,outputRoot,'Visible','off');
% Later, replot without rerunning inference:
% plot_kw51_modal_only(out,'Visible','on');
% Or pass a saved modal_... folder directly to plot_kw51_modal_only.
