function test_kw51_probabilistic(varargin)
% Independent numerical fixtures; run without original ZIPs.
% test_kw51_probabilistic('RunFit',true) additionally repeats one kinematics fit.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'functions'));
if exist('OCTAVE_VERSION','builtin')&&exist('butter','file')~=2,pkg('load','signal');end
ip=inputParser;addParameter(ip,'RunFit',false);parse(ip,varargin{:});
c=jsondecode(fileread(fullfile(root,'config.json')));kwpe_check_settings(c);
here=fileparts(mfilename('fullpath'));
F=load(fullfile(here,'core_fixture.mat'));
rng(19,'twister');c.posterior_draw_count=10000;
P=kwpe_condition(F.nom,F.delta,F.y,1:3,c,F.time_s);
assert(max(abs(P.weights(:)-F.expected_weights(:)))<1e-10,'Mixture-weight parity failed.');
assert(max(abs(P.mean(:)-F.expected_mean(:)))<1e-10,'Conditional-mean parity failed.');
assert(max(abs(P.structural_sd(:)-F.expected_sd(:)))<1e-10,'Conditional-SD parity failed.');
assert(abs(P.posterior_ess-F.expected_ess)<1e-10,'ESS parity failed.');
sampleVar=reshape(var(P.trajectories,0,1),size(P.mean));exactVar=P.structural_sd.^2;
relativeError=mean(abs(sampleVar(:)-exactVar(:))./max(exactVar(:),1e-12));
assert(relativeError<.05,'Conditional path variance mismatch.');
fprintf('PASS: independent Bridge West conditioning parity and path covariance (%.3f mean relative variance error).\n',relativeError);

F=load(fullfile(here,'fe_fixture.mat'));B=load(fullfile(root,'fe_models','model_15.mat'));
assert(size(B.G,2)==2,'Bundled FE bank is not the blueprint two-support bank.');
H=kw51_basis(B,F.E,F.E.time_s,F.F.kinematics);
assert(max(abs(H(:)-F.expected_basis(:)))<1e-8,'FE basis parity failed.');
c.candidate_count=12;c.paired_count=8;
[nom,delta]=kwpe_simulate(B,F.E,F.F,c,F.draws);
assert(max(abs(nom(:)-F.expected_nominal(:)))<1e-8,'Candidate FE response parity failed.');
assert(max(abs(delta(:)-F.expected_delta(:)))<1e-8,'Paired-FE discrepancy parity failed.');
fprintf('PASS: moving-load FE basis, candidate responses and paired discrepancy.\n');

E=F.E;kw51_check_input(E);
for name=[c.protected_names(:);c.same_station_exclusions(:);{'truth'}].'
    bad=E;bad.(name{1})=zeros(size(E.time_s));rejected=false;
    try,kw51_check_input(bad);catch,rejected=true;end
    assert(rejected,'Protected or undeclared input was accepted.');
end
fprintf('PASS: protected and same-station input fields rejected.\n');

K=load(fullfile(here,'filter_fixture.mat'));Y=kw51_filter(K.input,K.cutoff,K.fs);
assert(max(abs(Y(:)-K.expected(:)))<1e-9,'Filter parity failed.');
psi=kwpe_smooth_basis((0:.02:34.98).',.5);
assert(max(abs(sum(psi.^2,2)-1))<1e-12,'Smooth-process prior variance is not one.');
assert(kwpe_seed('traindata_20191210_170156',20260922)==F.expected_event_seed, ...
    'Event seed differs from the independent SHA-256 fixture.');
fprintf('PASS: filter, smooth-process variance and repeatable event seeding.\n');

if ip.Results.RunFit
    fitted=kw51_fit(B,F.E_full,1:3,8);
    assert(numel(fitted.kinematics)==4&&all(isfinite(fitted.kinematics)), ...
        'Fresh blueprint-BC kinematics fit returned invalid values.');
    assert(all(isfinite(fitted.loads_N))&&isfinite(fitted.objective), ...
        'Fresh blueprint-BC load/kinematics fit returned invalid values.');
    fprintf('PASS: fresh reference-only blueprint-BC kinematics fit smoke test.\n');
end
fprintf('All requested MATLAB self-tests passed. These tests do not fit protected outputs.\n');
end
