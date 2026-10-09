function out=run_kw51_modal_only(originalRoot,dataSource,outputRoot,varargin)
% Raw ZIPs -> modal-only predictions -> freeze -> withheld validation -> plots.
root=kwmo_setup(originalRoot);here=fileparts(mfilename('fullpath'));
ip=inputParser;addParameter(ip,'Events','all');addParameter(ip,'MaxEvents',Inf);
addParameter(ip,'Validate',true);addParameter(ip,'Plot',true);addParameter(ip,'Visible','off');
addParameter(ip,'RebuildModalBank',false);addParameter(ip,'CandidateCount',900);
addParameter(ip,'PairedCount',48);addParameter(ip,'PosteriorDrawCount',256);
addParameter(ip,'Seed',20260922);addParameter(ip,'PointsPerCycle',40);
addParameter(ip,'CheckEveryEvent',false);parse(ip,varargin{:});opt=ip.Results;
assert(exist('butter','file')==2,'Signal Processing Toolbox is required.');
assert(isscalar(opt.MaxEvents)&&opt.MaxEvents>=1&&(isinf(opt.MaxEvents)||opt.MaxEvents==floor(opt.MaxEvents)));
c=jsondecode(fileread(fullfile(here,'modal_only_config.json')));
c.seed=opt.Seed;c.candidate_count=opt.CandidateCount;c.paired_count=opt.PairedCount;
c.posterior_draw_count=opt.PosteriorDrawCount;c.points_per_cycle=opt.PointsPerCycle;
assert(c.candidate_count>=3&&c.candidate_count==floor(c.candidate_count));
assert(c.paired_count>=4&&c.paired_count==floor(c.paired_count));
assert(c.posterior_draw_count>=16&&c.posterior_draw_count==floor(c.posterior_draw_count));
assert(isfinite(c.seed)&&c.seed>=0&&c.seed==floor(c.seed));
assert(c.points_per_cycle>=20&&isfinite(c.points_per_cycle));
assert(c.baseline_stiffness_scale==1&&c.kinematic_support_force_weight==0);
assert(~isfield(c,'stiffness_states')&&~isfield(c,'pair_nominal_stiffness'));
if ~isfolder(outputRoot),mkdir(outputRoot);end
[ok,a]=fileattrib(outputRoot);assert(ok);outputRoot=a.Name;
base=fullfile(outputRoot,['modal_' datestr(now,'yyyymmdd_HHMMSS')]);run=base;i=0;
while isfolder(run),i=i+1;run=sprintf('%s_%03d',base,i);end
mkdir(run);for d={'predictions','draws','decomposition','fe_bank'},mkdir(fullfile(run,d{1}));end
bankFile=fullfile(run,'fe_bank','modal_bank.mat');
if opt.RebuildModalBank
 B=kwmo_build_modal_bank(root,bankFile,50);
else
 % The cached bank belongs to these exact FE source files.
 hashes=jsondecode(fileread(fullfile(here,'verified','fe_source_hashes.json')));
 for j=1:numel(hashes)
  path=fullfile(root,hashes(j).path);
  assert(isfile(path)&&strcmp(kw51_sha256(path),hashes(j).sha256), ...
   'Cached FE bank source mismatch: %s. Use RebuildModalBank=true to audit a rebuild.',path);
 end
 copyfile(fullfile(here,'verified','modal_bank.mat'),bankFile);B=load(bankFile);
end
baseB=load(fullfile(root,'fe_models','model_15.mat'));
assert(isequaln(B.H0,baseB.H0)&&isequaln(B.H2,baseB.H2)&&isequaln(B.G,baseB.G), ...
 'Modal bank differs from the supplied static blueprint bank.');
assert(isequal(B.labels(:),[c.reference_names(:);c.protected_names(:)]),'FE output order mismatch.');
staticBank=fullfile(run,'fe_bank','model_15.mat');kwpe_save_mat(staticBank,baseB);
csvFile=fullfile(run,'modal_samples.csv');copyfile(fullfile(here,'verified','modal_samples.csv'),csvFile);
S=kwmo_read_modal_samples(csvFile,B);kwpe_save_mat(fullfile(run,'modal_sample_metadata.mat'),S);
c.boundary_description=B.boundary_description;c.modal_modes=S.fe_modes;
c.modal_sample_count=numel(S.sample_ids);c.modal_source_sha256=S.source_sha256;
kwpe_save_mat(fullfile(run,'settings.mat'),c);
fid=fopen(fullfile(run,'settings.json'),'w');assert(fid>=0);fprintf(fid,'%s',jsonencode(c));fclose(fid);
fprintf('KW51 modal-only discrepancy: fixed unscaled blueprint FE, %d modal samples.\n%s\n',c.modal_sample_count,run);
catalog=kw51_catalog(dataSource,fullfile(outputRoot,'archive_cache'));
inputs=kw51_prepare(catalog,fullfile(run,'inputs'),opt.Events,opt.MaxEvents);
kwpe_save_mat(fullfile(run,'input_catalog.mat'),struct('inputs',inputs,'catalog',catalog));
rows={};
for k=1:numel(inputs)
 E=load(inputs(k).input_file);kw51_check_input(E);
 assert(max(S.sample_sdn)<E.strain_start_sdn,'Modal samples must precede every predicted passage.');
 timer=tic;fprintf('%d/%d: %s\n',k,numel(inputs),E.event_id);
 % Zero support-force weighting: reference-only kinematic fit, fixed FE.
 F=kw51_fit(B,E,1:3,0);
 [P,draws,D]=kwmo_predict(B,E,F,c,S,k==1||opt.CheckEveryEvent);
 if any(strcmp(c.development_events,E.event_id)),P.split='development';else,P.split='evaluation';end
 P.source_member=E.source_member;P.strain_start_sdn=E.strain_start_sdn;
 P.raw_archive=E.raw_archive;P.input_sha256=kw51_sha256(inputs(k).input_file);P.fit=F;P.seconds=toc(timer);
 kwpe_save_mat(fullfile(run,'predictions',[E.event_id '_prediction.mat']),P);
 kwpe_save_mat(fullfile(run,'draws',[E.event_id '_ensemble_draws.mat']),draws);
 kwpe_save_mat(fullfile(run,'decomposition',[E.event_id '_decomposition.mat']),D);
 rows(end+1,:)={E.event_id,P.split,E.track,F.kinematics(1),F.kinematics(2), ...
  numel(E.rail_offsets_s),P.posterior_ess,P.reference_nrmse(1),P.reference_nrmse(2),P.reference_nrmse(3),P.seconds}; %#ok<AGROW>
 fprintf('  ESS %.2f/%d; reference NRMSE %.3f %.3f %.3f; %.1fs\n',P.posterior_ess,c.candidate_count,P.reference_nrmse,P.seconds);
end
kw51_write_csv(fullfile(run,'reference_diagnostics.csv'), ...
 {'event','split','track','speed_m_s','entry_s','rail_peaks','ESS','NRMSE_ref1','NRMSE_ref2','NRMSE_ref3','seconds'},rows);
% Freeze inference and raw inputs before reading protected measurements.
kwpe_freeze(run,inputs,staticBank,root);
Fz=load(fullfile(run,'prediction_freeze.mat'));
extra=[kw51_files(fullfile(run,'decomposition'),'\.mat$'); ...
 {bankFile;csvFile;fullfile(run,'modal_sample_metadata.mat');fullfile(here,'modal_only_config.json')}];
code={'run_kw51_modal_only','kwmo_setup','kwmo_predict','kwmo_nominal','kwmo_modal_delta', ...
 'kwmo_newmark','kwmo_read_modal_samples','kwmo_decompose','kwmo_build_modal_bank','kwmo_solve_modes'};
for j=1:numel(code),extra{end+1,1}=fullfile(here,[code{j} '.m']);end
for j=1:numel(extra),Fz.files{end+1,1}=extra{j};Fz.hashes{end+1,1}=kw51_sha256(extra{j});end
kwpe_save_mat(fullfile(run,'prediction_freeze.mat'),Fz);
kw51_write_csv(fullfile(run,'prediction_freeze.csv'),{'file','SHA256'},[Fz.files,Fz.hashes]);
out=struct('run_dir',run,'output_root',outputRoot,'original_root',root, ...
 'events',{{inputs.event_id}},'settings',c,'fe_bank_file',bankFile,'validated',false);
kwpe_save_mat(fullfile(run,'run_output.mat'),struct('out',out));
if opt.Validate
 validate_kw51_probabilistic_matlab(run);out.validated=true;
 kwpe_save_mat(fullfile(run,'run_output.mat'),struct('out',out));
 quantify_kw51_modal_discrepancy(out);
end
if opt.Plot&&out.validated,plot_kw51_modal_only(out,'Visible',opt.Visible);end
fprintf('Completed modal-only run: %s\n',run);
end
