function outputRoot=run_bridgewest_stabil_shaker_informed(dataRoot,outputRoot,varargin)
%RUN_BRIDGEWEST_STABIL_SHAKER_INFORMED Reference-only 36-passage prediction.
%
% Scientific protocol
% -------------------
% For each requested train passage, ONLY measured A031 and B053 are read as
% structural inference channels. They are used to infer the eight bounded
% bogie resultants and condition the frozen 84-state shaker-informed Stabil
% discrepancy. The model then predicts A058, A070, B031 and B081.
%
% Protected target measurements are NOT read/preprocessed by this function.
% They are opened only by score_bridgewest_stabil_results after prediction.
%
% Default event set: 12 passages from 2022 + 24 passages from 2025 = 36.
%
% Example:
%   run_bridgewest_stabil_shaker_informed('D:\BridgeWestData','D:\BW_Stabil_v3');
root=fileparts(mfilename('fullpath'));pkg=fileparts(fileparts(root));
addpath(root,fullfile(root,'raw_functions'),fullfile(pkg,'matlab'), ...
    fullfile(pkg,'shaker_informed','matlab'),fullfile(pkg,'dependencies','Stabil'));
assert(nargin>=1&&~isempty(dataRoot),'Pass the folder containing LT-W ZIP archives and New folder (2).zip.');
if nargin<2||isempty(outputRoot),outputRoot=fullfile(root,'results_all_passages');end
ip=inputParser;
addParameter(ip,'Campaign','all');
addParameter(ip,'Events',{});
addParameter(ip,'ContinueOnError',false,@islogical);
parse(ip,varargin{:});opt=ip.Results;
assert(exist('butter','file')==2&&exist('filtfilt','file')==2,'Raw preprocessing requires Signal Processing Toolbox.');
if ~exist(outputRoot,'dir'),mkdir(outputRoot);end
for d={'processed_references','permitted_inputs','predictions','diagnostics','fe_model','tables'}
    q=fullfile(outputRoot,d{1});if ~exist(q,'dir'),mkdir(q);end
end
pre=jsondecode(fileread(fullfile(pkg,'production','config','preprocessing_config.json')));
c=jsondecode(fileread(fullfile(pkg,'production','config','prediction_config.json')));
evset=jsondecode(fileread(fullfile(pkg,'production','config','event_sets.json')));
sc=jsondecode(fileread(fullfile(pkg,'shaker_informed','discrepancy_config.json')));
fields={'cov_ridge','reference_noise_sigma_microstrain','target_noise_sigma_microstrain','displacement_noise_sigma_mm','interval_z'};
for j=1:numel(fields),sc.(fields{j})=c.(fields{j});end
A=bw_prediction_zip_index(dataRoot,fullfile(outputRoot,'archive_cache'));
assert(A.counts.ltw2022>=26&&A.counts.ltw2025>=50,'Unexpected LT-W archive inventory.');

% Detailed Stabil FE and shaker ensemble are fixed before traffic inference.
fe=run_bridgewest_lorenzen(fullfile(pkg,'config','model.json'),fullfile(outputRoot,'fe_model'));
% Exact 41-point BA/BB displacement coordinates from the detailed observation map.
if isfield(fe,'observation_maps') && ~isempty(fe.observation_maps)
    maps=fe.observation_maps;
else
    maps=bw_lorenzen_observation_maps(fe);
end
dispX=double(maps.x_dataset_m(:)).';
assert(numel(dispX)==41,'Expected 41 detailed-model displacement stations per girder.');
ens=bw_shaker_ensemble(fullfile(pkg,'shaker_informed'));
weights=double(ens.states.weight);weights=weights/sum(weights);
events=resolveEvents(opt,evset);
refsMeasured=cellstr(string(c.measured_reference_channels));
targetsMeasured=cellstr(string(c.measured_target_channels));
status=struct('event',{},'status',{},'reason',{});
manifest=struct('event',{},'file',{},'sha256',{},'speed_m_s',{},'load_condition_number',{});
settings=struct( ...
    'package_version','BridgeWest_Stabil_ShakerInformed_Discrepancy_v3_1', ...
    'protocol','A031/B053-only passage inference; A058/A070/B031/B081 prediction; target measurements assessment-only', ...
    'events',{events},'n_requested_passages',numel(events), ...
    'references',{refsMeasured},'protected_targets',{targetsMeasured}, ...
    'prediction_config',c,'preprocessing_config',pre,'shaker_summary',ens.summary, ...
    'model_geometry_verified',fe.c.geometry_verified,'table_5_2_imported',false, ...
    'target_histories_read_during_inference',false, ...
    'additional_target_derived_discrepancy',false, ...
    'displacement_prediction_active',true, ...
    'displacement_station_count_per_girder',numel(dispX));
bw_write_json(fullfile(outputRoot,'run_settings.json'),settings);

for ie=1:numel(events)
    eid=events{ie};fprintf('\n[%d/%d] %s\n',ie,numel(events),eid);
    try
        % This preprocessing call receives ONLY A031/B053.
        z=bw_prepare_raw_event(A,eid,pre,outputRoot,refsMeasured);
        [cols,reason]=fixedMeasuredReferences(z,refsMeasured);
        if ~isempty(reason)
            status(end+1)=struct('event',eid,'status','skipped','reason',reason); %#ok<AGROW>
            continue;
        end
        allowed=struct('event',eid,'time_s',z.time(:),'speed_m_s',z.speed,'entry_s',z.entry, ...
            'reference_names',{refsMeasured},'reference_strain_microstrain',z.Y(:,cols), ...
            'source',z.source,'contains_protected_target_measurements',false);
        save(fullfile(outputRoot,'permitted_inputs',[eid '.mat']),'-struct','allowed','-v7');

        % Only the permitted object is retained beyond this point.
        t=allowed.time_s;v=allowed.speed_m_s;entry=allowed.entry_s;yref=allowed.reference_strain_microstrain;
        clear z allowed

        [B0,meta0]=bw_stabil_basis_shaker(fe,v,ens.baseline,t,sc,1);
        pair=bw_stabil_paired_discrepancy(fe,v,t,[],sc,1,'');
        refs=double(c.reference_model_indices(:)).';
        [p,dg]=bw_predict_continuous_shaker(B0,pair.delta_basis,weights,yref,refs,c);

        pred=struct('event',eid,'time_s',t,'speed_m_s',v,'entry_s',entry, ...
            'model_channel_order',{cellstr(string(c.model_channel_order))}, ...
            'reference_model_indices',refs,'target_model_indices',double(c.target_model_indices(:)).', ...
            'target_labels',{cellstr(string(c.target_labels))}, ...
            'strain_mean_microstrain',p.strain_mean,'strain_sd_microstrain',p.strain_sd, ...
            'strain_lo_microstrain',p.strain_lo,'strain_hi_microstrain',p.strain_hi, ...
            'disp_BA_mean_mm',p.disp_BA_mean,'disp_BA_sd_mm',p.disp_BA_sd, ...
            'disp_BA_lo_mm',p.disp_BA_lo,'disp_BA_hi_mm',p.disp_BA_hi, ...
            'disp_BB_mean_mm',p.disp_BB_mean,'disp_BB_sd_mm',p.disp_BB_sd, ...
            'disp_BB_lo_mm',p.disp_BB_lo,'disp_BB_hi_mm',p.disp_BB_hi, ...
            'displacement_x_dataset_m',dispX, ...
            'load_mean_kN',p.load_mean_kN,'load_sd_kN',p.load_sd_kN, ...
            'load_covariance_kN2',p.load_covariance_kN2,'load_correlation',p.load_correlation, ...
            'load_condition_number',p.load_condition_number, ...
            'shaker_effective_state_count',p.shaker_effective_state_count, ...
            'fixed_point_trace',{p.iteration_trace},'nominal_modal_state',meta0.modal_state, ...
            'protected_targets_used_in_inference',false, ...
            'target_measurements_loaded_during_inference',false, ...
            'run40_status',char(evset.run40_status));
        file=fullfile(outputRoot,'predictions',[eid '.mat']);save(file,'-struct','pred','-v7.3');
        dsave=struct('event',eid,'load_posterior',dg.load_posterior, ...
            'fixed_point_load_kN',dg.fixed_point_load_kN,'trace',{dg.trace});
        save(fullfile(outputRoot,'diagnostics',[eid '_diagnostics.mat']),'-struct','dsave','-v7');
        manifest(end+1)=struct('event',eid,'file',file,'sha256',bw_sha256(file), ...
            'speed_m_s',v,'load_condition_number',p.load_condition_number); %#ok<AGROW>
        status(end+1)=struct('event',eid,'status','predicted', ...
            'reason','A031/B053 only; A058/A070/B031/B081 not read during inference'); %#ok<AGROW>
        fprintf('  p(kN) = %s\n',sprintf('%.1f ',p.load_mean_kN));
    catch ME
        status(end+1)=struct('event',eid,'status','failed','reason',ME.message); %#ok<AGROW>
        fprintf(2,'  FAILED: %s\n',ME.message);
        if ~opt.ContinueOnError,rethrow(ME);end
    end
    bw_write_json(fullfile(outputRoot,'event_status.json'),status);
end
bw_write_json(fullfile(outputRoot,'prediction_manifest.json'),manifest);

nPred=sum(strcmp({status.status},'predicted'));
nSkip=sum(strcmp({status.status},'skipped'));
nFail=sum(strcmp({status.status},'failed'));
stageSummary=struct('requested_passages',numel(events),'predicted_passages',nPred, ...
    'skipped_passages',nSkip,'failed_passages',nFail,'predicted_target_histories',4*nPred, ...
    'predicted_displacement_histories',2*nPred, ...
    'displacement_station_count_per_girder',numel(dispX), ...
    'measured_structural_inference_channels',{refsMeasured}, ...
    'withheld_target_channels',{targetsMeasured}, ...
    'target_measurements_read_during_inference',false);
bw_write_json(fullfile(outputRoot,'prediction_stage_summary.json'),stageSummary);

isDefaultAll=isempty(opt.Events)&&lower(string(opt.Campaign))=="all";
if isDefaultAll&&~opt.ContinueOnError
    assert(numel(events)==36,'Default all-passage event set must contain 36 passages.');
    assert(nPred==36&&nSkip==0&&nFail==0, ...
        'Reference-only production run did not produce all 36 passage predictions. See event_status.json.');
end
fprintf('\nReference-only prediction stage complete: %d/%d passages predicted.\n',nPred,numel(events));
fprintf('Protected target measurements have NOT been read.\n');
fprintf('Run score_bridgewest_stabil_results(dataRoot,''%s'') as the separate assessment step.\n',outputRoot);
end

function events=resolveEvents(opt,evset)
if ~isempty(opt.Events)
    if ischar(opt.Events)||isstring(opt.Events),events=cellstr(string(opt.Events));else,events=cellstr(string(opt.Events(:)));end
    return;
end
c=lower(string(opt.Campaign));events={};
if c=="all"||c=="2022"
    for r=double(evset.runs_2022(:)).',events{end+1}=sprintf('2022_run%02d',r);end %#ok<AGROW>
end
if c=="all"||c=="2025"
    for r=double(evset.runs_2025(:)).',events{end+1}=sprintf('2025_run%02d',r);end %#ok<AGROW>
end
assert(~isempty(events),'Campaign must be all, 2022 or 2025, or supply Events.');
end

function [cols,reason]=fixedMeasuredReferences(z,expected)
cols=zeros(1,numel(expected));reason='';chs=cellstr(string(z.channels));flags=cellstr(string(z.flags));
for j=1:numel(expected)
    ii=find(strcmp(chs,expected{j}));
    if numel(ii)~=1,reason=sprintf('Required reference %s missing or duplicated.',expected{j});cols=[];return;end
    if ~strcmp(flags{ii},'ok')||~all(isfinite(z.Y(:,ii))),reason=sprintf('Required reference %s invalid.',expected{j});cols=[];return;end
    cols(j)=ii;
end
end
