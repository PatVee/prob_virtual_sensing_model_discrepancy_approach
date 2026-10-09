function outDir=run_bridgewest_known_load_diagnostic(outputRoot,varargin)
%RUN_BRIDGEWEST_KNOWN_LOAD_DIAGNOSTIC Fixed-FE, static-weighing comparison.
% Optional diagnostic, NOT a replacement for the reference-conditioned solver.
% For each selected passage, regenerates the nominal unit-bogie FE basis and
% applies (a) saved inferred bogie means and (b) original static weighed bogies.
% Uses no target strain or PLS data. Both strain and displacement are propagated.
% No scalar correction is fitted to measured displacement; mass scale stays 1.
%
%   run_bridgewest_known_load_diagnostic(outputRoot,'Event','2025_run22');
%   run_bridgewest_known_load_diagnostic(outputRoot,'Event','all');
%
% These are deterministic NOMINAL-FE responses, without reference-residual
% conditioning. They must not be confused with the saved conditional mean.
% Static weighed loads do not constitute measured dynamic contact histories.
ip=inputParser;addParameter(ip,'Event',[]);addParameter(ip,'Campaign','');addParameter(ip,'Runs',[]);
addParameter(ip,'OutputDir','');addParameter(ip,'WeightDir','');parse(ip,varargin{:});o=ip.Results;
here=fileparts(mfilename('fullpath'));pkg=fileparts(fileparts(here));
addpath(here,fullfile(pkg,'matlab'),fullfile(pkg,'shaker_informed','matlab'),fullfile(pkg,'dependencies','Stabil'));
events=bw_select_saved_events(outputRoot,'Event',o.Event,'Campaign',o.Campaign,'Runs',o.Runs);
if isempty(o.OutputDir),outDir=fullfile(outputRoot,'known_load_diagnostic');else,outDir=o.OutputDir;end
assert(~strcmp(char(outDir),fullfile(char(outputRoot),'predictions')),'BridgeWest:ProtectedPath','Do not overwrite the production predictions directory.');
if ~isfolder(outDir),mkdir(outDir);end
S=load(fullfile(outputRoot,'fe_model','bridgewest_lorenzen_fe.mat'),'fe');fe=S.fe;
c=jsondecode(fileread(fullfile(pkg,'shaker_informed','discrepancy_config.json')));
for i=1:numel(events)
    eid=events{i};P=load(fullfile(outputRoot,'predictions',[eid '.mat']));
    W=bw_weighed_bogie_loads(str2double(eid(1:4)),o.WeightDir);
    assert(isfield(P,'nominal_modal_state'),'BridgeWest:MissingState','Saved nominal shaker state is missing for %s.',eid);
    % v3/v3.1 simulated the requested DT passage set with direction=+1.
    % Reuse this convention; do not infer a direction from target waveforms.
    [B,meta]=bw_stabil_basis_shaker(fe,P.speed_m_s,P.nominal_modal_state,P.time_s,c,1);
    inferred=double(P.load_mean_kN(:))*1000;weighed=W.bogie_load_kN(:)*1000;
    O=struct('event',eid,'time_s',P.time_s,'speed_m_s',P.speed_m_s, ...
        'model_channel_order',{P.model_channel_order},'displacement_x_dataset_m',P.displacement_x_dataset_m, ...
        'inferred_bogie_kN',inferred.'/1000,'weighed_bogie_kN',weighed.'/1000, ...
        'bridge_mass_scale',1,'target_measurements_used',false,'PLS_used',false, ...
        'reference_residual_conditioning_applied',false, ...
        'description','Nominal FE diagnostic only; static weighed inputs; no interval or accuracy claim', ...
        'weighing_source',W.source_file,'basis_meta',meta);
    fields={'strain_microstrain','disp_BA_mm','disp_BB_mm'};
    for j=1:3
        nt=size(B{j},1);nc=size(B{j},2);bb=reshape(B{j},nt*nc,8);
        O.(['inferred_load_FE_' fields{j}])=reshape(bb*inferred,nt,nc);
        O.(['weighed_load_FE_' fields{j}])=reshape(bb*weighed,nt,nc);
    end
    save(fullfile(outDir,[eid '_known_load.mat']),'-struct','O','-v7.3');
    fprintf('[%d/%d] %s: inferred %.2f kN; weighed %.2f kN. No mass adjustment.\n', ...
        i,numel(events),eid,sum(inferred)/1000,sum(weighed)/1000);
end
end
