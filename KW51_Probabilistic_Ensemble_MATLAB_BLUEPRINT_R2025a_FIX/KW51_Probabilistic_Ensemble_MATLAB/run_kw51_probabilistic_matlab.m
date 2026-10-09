function out=run_kw51_probabilistic_matlab(dataSource,outputRoot,varargin)
%RUN_KW51_PROBABILISTIC_MATLAB Full blueprint-BC probabilistic workflow.
%
% Original ZIP/raw MAT files -> fresh passage inference -> validation.
% The nominal FE bank uses the official KW51 blueprint boundary condition:
% Brussels longitudinally anchored, Leuven longitudinal DOFs free, kxx=0.
% The paired support discrepancy acts only at the two Leuven longitudinal
% bearing coordinates. No independent common-translation prior is used.

root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'functions'));
if nargin<1||isempty(dataSource),dataSource=fullfile(root,'data');end
if nargin<2||isempty(outputRoot),outputRoot=fullfile(root,'KW51ProbabilisticResults');end
if exist('OCTAVE_VERSION','builtin')&&exist('butter','file')~=2
    try,pkg('load','signal');catch,error('The signal package is required.');end
end
assert(exist('butter','file')==2,'Signal Processing Toolbox (butter) is required.');

ip=inputParser;
addParameter(ip,'Events','all');
addParameter(ip,'MaxEvents',Inf);
addParameter(ip,'Validate',true);
addParameter(ip,'Plot',true);
addParameter(ip,'RebuildFE',false);
addParameter(ip,'Seed',20260922);
addParameter(ip,'CandidateCount',900);
addParameter(ip,'PairedCount',48);
addParameter(ip,'PosteriorDrawCount',256);
% Retained only for API compatibility. Blueprint BC fixes this to zero.
addParameter(ip,'TranslationSDmm',0);
parse(ip,varargin{:});
opt=ip.Results;

assert(isscalar(opt.MaxEvents)&&opt.MaxEvents>0&& ...
    (isinf(opt.MaxEvents)||opt.MaxEvents==floor(opt.MaxEvents)), ...
    'MaxEvents must be a positive integer or Inf.');
assert(isscalar(opt.TranslationSDmm)&&abs(opt.TranslationSDmm)<1e-12, ...
    ['Blueprint BC does not use an independent common-translation prior. ' ...
     'TranslationSDmm must be 0.']);

c=jsondecode(fileread(fullfile(root,'config.json')));
c.seed=opt.Seed;
c.candidate_count=opt.CandidateCount;
c.paired_count=opt.PairedCount;
c.posterior_draw_count=opt.PosteriorDrawCount;
c.translation_sd_mm=0;
c.translation_sensitivity_mm=0;
kwpe_check_settings(c);

if ~exist(outputRoot,'dir'),mkdir(outputRoot);end
[ok,info]=fileattrib(outputRoot);assert(ok);outputRoot=info.Name;
base=fullfile(outputRoot,['run_' datestr(now,'yyyymmdd_HHMMSS')]);
runFolder=base;k=0;
while exist(runFolder,'dir'),k=k+1;runFolder=sprintf('%s_%03d',base,k);end
mkdir(runFolder);
mkdir(fullfile(runFolder,'predictions'));
mkdir(fullfile(runFolder,'draws'));
mkdir(fullfile(runFolder,'fe_bank'));

diary(fullfile(runFolder,'execution_log.txt'));
diaryCleanup=onCleanup(@()diary('off')); %#ok<NASGU>
fprintf('KW51 probabilistic ensemble - BLUEPRINT BC - fresh native MATLAB execution\n%s\n',runFolder);
fprintf('Protected outputs: sgBD1718C, sgBD2324C, dABLA, dABLC. Validation only.\n');
fprintf('No independent common longitudinal translation prior is used.\n');

% Every run owns an exact local FE bank. Plotting/post-processing therefore
% cannot accidentally replay draws with a different boundary condition.
runBankFolder=fullfile(runFolder,'fe_bank');
runBankFile=fullfile(runBankFolder,'model_15.mat');
if opt.RebuildFE
    build_blueprint_bank(runBankFolder);
else
    bundledBank=fullfile(root,'fe_models','model_15.mat');
    assert(exist(bundledBank,'file')==2, ...
        'Bundled blueprint FE bank is missing. Use RebuildFE=true.');
    copyfile(bundledBank,runBankFile,'f');
end
assert(exist(runBankFile,'file')==2,'Run-local FE response bank was not created.');
B=load(runBankFile);
validate_blueprint_bank(B,c);

c.boundary_description=B.boundary_description;
c.support_coordinate_count=size(B.G,2);
kwpe_save_mat(fullfile(runFolder,'settings.mat'),c);
fid=fopen(fullfile(runFolder,'settings.json'),'w');assert(fid>=0);
fprintf(fid,'%s',jsonencode(c));fclose(fid);

fprintf('FE bank: %s\n',runBankFile);
fprintf('BC: %s\n',B.boundary_description);
fprintf('Support discrepancy coordinates: %d\n',size(B.G,2));

catalog=kw51_catalog(dataSource,fullfile(outputRoot,'archive_cache'));
inputs=kw51_prepare(catalog,fullfile(runFolder,'inputs'),opt.Events,opt.MaxEvents);
kwpe_save_mat(fullfile(runFolder,'input_catalog.mat'),struct('inputs',inputs,'catalog',catalog));

diagnostics={};
for ie=1:numel(inputs)
    E=load(inputs(ie).input_file);kw51_check_input(E);timer=tic;
    fprintf('%d/%d: %s, fresh kinematic fit\n',ie,numel(inputs),E.event_id);
    % Auxiliary load fit estimates kinematics only. Its fitted point loads
    % are not substituted for the independent probabilistic load ensemble.
    F=kw51_fit(B,E,1:3,c.support_force_sd_N/1e5);
    [P,draws]=kwpe_predict(B,E,F,c);
    if any(strcmp(c.development_events,E.event_id)),P.split='development';else,P.split='evaluation';end
    P.source_member=E.source_member;P.strain_start_sdn=E.strain_start_sdn;
    P.raw_archive=E.raw_archive;P.input_sha256=kw51_sha256(inputs(ie).input_file);
    P.fit=F;P.seconds=toc(timer);
    kwpe_save_mat(fullfile(runFolder,'predictions',[E.event_id '_prediction.mat']),P);
    kwpe_save_mat(fullfile(runFolder,'draws',[E.event_id '_ensemble_draws.mat']),draws);
    diagnostics(end+1,:)={E.event_id,P.split,E.track,F.kinematics(1),F.kinematics(2), ...
        numel(E.rail_offsets_s),P.posterior_ess,P.reference_nrmse(1),P.reference_nrmse(2), ...
        P.reference_nrmse(3),P.seconds}; %#ok<AGROW>
    fprintf('  ESS %.2f / %d; reference NRMSE [%.3f %.3f %.3f]; %.1f s\n', ...
        P.posterior_ess,c.candidate_count,P.reference_nrmse,P.seconds);
end
kw51_write_csv(fullfile(runFolder,'reference_diagnostics.csv'), ...
    {'event','split','track','speed_m_s','entry_s','rail_peaks','ESS', ...
    'NRMSE_ref1','NRMSE_ref2','NRMSE_ref3','seconds'},diagnostics);

kwpe_freeze(runFolder,inputs,runBankFile,root);
out=struct('run_dir',runFolder,'output_root',outputRoot,'events',{{inputs.event_id}}, ...
    'settings',c,'fe_bank_file',runBankFile,'validated',false);
fprintf('All predictions saved and frozen before target validation.\n');
if opt.Validate
    validate_kw51_probabilistic_matlab(runFolder);out.validated=true;
end
kwpe_save_mat(fullfile(runFolder,'run_output.mat'),struct('out',out));
if opt.Plot
    if out.validated
        plot_kw51_probabilistic_matlab(out);
    else
        warning('KW51:NoValidation','Plotting measured comparisons was skipped because Validate=false.');
    end
end
fprintf('Completed: %s\n',runFolder);
end

function validate_blueprint_bank(B,c)
required={'H0','H2','G','xgrid','labels','boundary_description'};
for k=1:numel(required)
    assert(isfield(B,required{k}),'FE bank is missing field %s.',required{k});
end
assert(isequal(B.labels(:),[c.reference_names(:);c.protected_names(:)]), ...
    'FE channel order mismatch.');
assert(isequal(size(B.G),[7 2]), ...
    ['This package requires the blueprint FE bank with exactly two Leuven ' ...
     'longitudinal support coordinates.']);
txt=lower(char(B.boundary_description));
assert(contains(txt,'blueprint')&&contains(txt,'leuven'), ...
    'The selected FE bank is not identified as the blueprint BC bank.');
if isfield(B,'common_translation_removed')
    assert(~logical(B.common_translation_removed), ...
        'A legacy zero-mean floating bank was supplied instead of the blueprint bank.');
end
end
