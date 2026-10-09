function ens = bw_shaker_ensemble(root)
%BW_SHAKER_ENSEMBLE Load the shaker-informed modal-state mixture.
% 2024 contributes genuine same-run bending/torsion pairs. 2025 contributes
% the product of the two empirical mode marginals within each excitation level
% because the modes were measured in separate run blocks; no rank-pairing is used.
if nargin<1 || isempty(root)
    root=fileparts(fileparts(mfilename('fullpath')));
end
res=fullfile(root,'results');
ens.states=readtable(fullfile(res,'shaker_weighted_states_all.csv'),'TextType','string');
ens.joint_2024=readtable(fullfile(res,'shaker_joint_states_2024.csv'),'TextType','string');
ens.product_2025=readtable(fullfile(res,'shaker_product_states_2025.csv'),'TextType','string');
ens.crosscheck_2025=readtable(fullfile(res,'shaker_crosscampaign_2025_unpaired.csv'),'TextType','string');
ens.summary=jsondecode(fileread(fullfile(res,'shaker_ensemble_summary.json')));
b=ens.summary.baseline;
ens.baseline=struct('f1_Hz',b.mode1_frequency_Hz,'f2_Hz',b.mode2_frequency_Hz, ...
    'zeta1',b.mode1_damping_ratio,'zeta2',b.mode2_damping_ratio,'label','pooled_5mps2_median');
assert(abs(sum(ens.states.weight)-1)<1e-10,'Shaker state weights do not sum to one.');
end
