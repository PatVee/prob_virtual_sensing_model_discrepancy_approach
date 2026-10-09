%RUN_SHAKER_INFORMED_DISCREPANCY_DEMO End-to-end mechanics/discrepancy demonstration.
% This does not use withheld train target strains or PLS displacement.
clear;clc;
root=fileparts(fileparts(mfilename('fullpath')));pkg=fileparts(root);
addpath(fullfile(pkg,'matlab'));addpath(fullfile(root,'matlab'));addpath(fullfile(pkg,'dependencies','Stabil'));
c=jsondecode(fileread(fullfile(root,'discrepancy_config.json')));
fe=run_bridgewest_lorenzen(fullfile(pkg,'config','model.json'),fullfile(pkg,'results','matlab'));
v=160/3.6;span=diff(fe.c.bearing_x_m);
tobs=(-c.pre_entry_s:1/50:(span+max(c.axle_offsets_m))/v+c.post_exit_s).';
% Use one common mean-load vector so the demo covariance reflects only the
% shaker-state dynamics. Production inference may re-project out.delta_basis
% with a passage-specific posterior mean or candidate load vector.
meanLoad=mean(double(c.default_bogie_load_range_N));loadsN=ones(1,8)*meanLoad;
out=bw_stabil_paired_discrepancy(fe,v,tobs,loadsN,c,1,fullfile(root,'results','demo_paired_discrepancy.mat'));
refs=double(c.reference_indices_matlab(:)).';
pc=struct('cov_ridge',c.cov_ridge,'reference_noise_sigma_microstrain',c.reference_noise_sigma_microstrain, ...
    'target_noise_sigma_microstrain',c.target_noise_sigma_microstrain,'displacement_noise_sigma_mm',c.displacement_noise_sigma_mm);
disc=bw_simulation_discrepancy_weighted(out.delta,out.weights,refs,pc);
save(fullfile(root,'results','demo_discrepancy_moments.mat'),'disc','-v7.3');
fprintf('Built %d weighted shaker states (effective count %.2f).\n',disc.paired_draw_count,disc.effective_state_count);
fprintf('Reference channels: A031 and B053. No withheld target was used.\n');
