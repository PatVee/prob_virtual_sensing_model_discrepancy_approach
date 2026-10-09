function build_blueprint_bank(outputFolder)
%BUILD_BLUEPRINT_BANK Build the KW51 quasi-static response bank using the
% official blueprint support conditions.
%
% Nominal support definition (official model, fixed=0, kxx=0):
%   Brussels N1: x,y,z restrained
%   Brussels N2: x,z restrained
%   Leuven   N3: y,z restrained; x free
%   Leuven   N4: z restrained; x,y free
%
% The two equivalent support-force discrepancy coordinates are the
% longitudinal Leuven bearing DOFs corresponding to dABLA and dABLC.
% No four-bearing floating release and no zero-mean rigid-translation gauge
% are introduced.

root=fileparts(mfilename('fullpath'));
fe=fullfile(root,'dependencies','official_fe');
addpath(genpath(fullfile(root,'dependencies','Stabil')));
addpath(fe);
if nargin<1||isempty(outputFolder),outputFolder=fullfile(root,'fe_models');end
if ~exist(outputFolder,'dir'),mkdir(outputFolder);end

old=pwd;cleanup=onCleanup(@()cd(old)); %#ok<NASGU>
cd(fe);

S=load('sensor_definition.mat','Ss','Sd','labels_s','labels_d');
labels={'sgBD1011A','sgBD1415A','sgBD2728A', ...
    'sgBD1718C','sgBD2324C','dABLA','dABLC'};
[ok,ix]=ismember(labels(1:5),S.labels_s);
assert(all(ok),'One or more required strain channels are absent from sensor_definition.mat.');
assert(isequal(S.labels_d(:),{'dABLA';'dABLC'}),'Unexpected displacement-channel order.');
H=[S.Ss(ix,:)*1e6;S.Sd*1e3]; % microstrain and mm

% Keep the existing ensemble nominal structural parameters; change only BC.
% model_kw51 inputs: Es_rel,hb_rel,kxx,Id_rel,Ad_rel,fixed,plot,addmass,diags
[K,~,T]=model_kw51(1,1,0,1,1,0,0,1,[]);

% Full-force -> selected-output static operator.
Z=K\(T.'*H.');
Q=Z.'*T.';

% Two physically retained longitudinal support discrepancy coordinates.
% Since d = Sd*U, virtual work gives F_full = Sd'*f_support.
G=full(Q*S.Sd.');
assert(isequal(size(G),[7 2]),'Blueprint support map must be 7 x 2.');

xgrid=(-5:.5:121.5).';
tracks={'A','B'};
H0=zeros(numel(xgrid),7,2);H2=H0;
for tr=1:2
    P0=xrail_2_P(tracks{tr},xgrid,1,zeros(size(xgrid)));
    P30=xrail_2_P(tracks{tr},xgrid,1,30*ones(size(xgrid)));
    P2=(P30-P0)/900;
    H0(:,:,tr)=(Q*P0).';
    H2(:,:,tr)=(Q*P2).';
end

parameter=[1 1 1];
compliance_m_per_N=G(6:7,:)/1000;
boundary_description=['KW51 blueprint BC: Brussels longitudinally anchored; ' ...
    'Leuven longitudinal DOFs free; kxx=0.'];
support_force_names={'Leuven_A_longitudinal','Leuven_C_longitudinal'};
common_translation_removed=false;
blueprint_fixed=0;
leuven_kxx_N_m=0;

save(fullfile(outputFolder,'model_15.mat'), ...
    'H0','H2','G','compliance_m_per_N','parameter','xgrid','labels', ...
    'boundary_description','support_force_names','common_translation_removed', ...
    'blueprint_fixed','leuven_kxx_N_m','-v7');

fprintf('\nSaved KW51 BLUEPRINT response bank.\n');
fprintf('  fixed                         = 0\n');
fprintf('  Leuven kxx                   = 0 N/m\n');
fprintf('  support discrepancy columns = %d\n',size(G,2));
fprintf('  supports                    = N3-x, N4-x (Leuven)\n\n');
end
