function test_kw51_modal_only(originalRoot)
% Tests specifically guard removal of stiffness scaling and missing plot data.
kwmo_setup(originalRoot);here=fileparts(mfilename('fullpath'));
c=jsondecode(fileread(fullfile(here,'modal_only_config.json')));
assert(~isfield(c,'stiffness_states')&&~isfield(c,'pair_nominal_stiffness'));
assert(c.baseline_stiffness_scale==1&&c.kinematic_support_force_weight==0);
B=struct('frequency_Hz',2,'xgrid',[-5;0;5;10;121.5], ...
 'modal_force0',ones(5,1,2),'modal_force2',zeros(5,1,2),'output_modes',ones(7,1));
E=struct('time_s',(0:.02:2).','rail_offsets_s',0,'rail_first_time_s',.5,'track','A');
fit=struct('kinematics',[10 .5 1 0]);
S=struct('fe_modes',1,'median_damping',.03,'frequency_ratio',ones(3,1), ...
 'damping_ratio',.03*ones(3,1),'source_sha256','synthetic');
rng(7);[d,draws]=kwmo_modal_delta(B,E,fit,c,ones(3,1),S);
assert(max(abs(d(:)))==0,'Identical modal parameters must give zero discrepancy.');
assert(draws.nominal_frequency_Hz==2,'The old sqrt(0.9) scaling has reappeared.');
% Analytical integration check.
h=.001;t=(0:h:2).';z=.03;w=4*pi;wd=w*sqrt(1-z*z);
q=kwmo_newmark(ones(size(t)),h,2,z);
exact=(1-exp(-z*w*t).*(cos(wd*t)+z/sqrt(1-z*z)*sin(wd*t)))/w^2;
assert(max(abs(q-exact))<2e-6);
% Missing channels are valid; matrices are not silently flattened.
V=struct('native_time',{{[]}},'native_truth',{{[]}},'native_prediction',{{[]}}, ...
 'native_lower95',{{[]}},'native_upper95',{{[]}},'native_score_mask',{{[]}});
A=kwmo_native_series(V,1);assert(~A.available);
V.native_time{1}=[0 1 2];V.native_truth{1}=[1;2;3];V.native_prediction{1}=[1 2 3];
A=kwmo_native_series(V,1);assert(A.available&&isequal(size(A.time),[3 1]));
V.native_truth{1}=ones(3,2);caught=false;
try,kwmo_native_series(V,1);catch,caught=true;end
assert(caught,'Matrix-valued measurements should be rejected.');
B=load(fullfile(here,'verified','modal_bank.mat'));
S=kwmo_read_modal_samples(fullfile(here,'verified','modal_samples.csv'),B);
assert(size(S.frequency_ratio,1)==916&&numel(S.fe_modes)==6);
assert(max(abs(median(S.frequency_ratio,1)-1))<1e-12);
fprintf('PASS: unscaled modal centre, zero perturbation, integrator, missing-channel handling and sample reader.\n');
end
