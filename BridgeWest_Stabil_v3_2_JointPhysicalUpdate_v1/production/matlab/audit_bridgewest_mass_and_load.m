function A=audit_bridgewest_mass_and_load(outputRoot,varargin)
%AUDIT_BRIDGEWEST_MASS_AND_LOAD Inspect saved FE mass and inferred load scale.
% Reads saved FE/predictions and original weighing CSVs. Does not read PLS,
% update the FE, change inferred loads, or modify any prediction file.
%   A=audit_bridgewest_mass_and_load(outputRoot);
ip=inputParser;addParameter(ip,'OutputDir','');addParameter(ip,'WeightDir','');parse(ip,varargin{:});o=ip.Results;
if isempty(o.OutputDir),o.OutputDir=fullfile(outputRoot,'audit_mass_load');end
if ~isfolder(o.OutputDir),mkdir(o.OutputDir);end
ff=fullfile(outputRoot,'fe_model','bridgewest_lorenzen_fe.mat');
assert(isfile(ff),'BridgeWest:MissingFE','Saved native FE model is missing: %s',ff);
S=load(ff,'fe');assert(isfield(S,'fe'),'BridgeWest:FE','Expected fe structure.');fe=S.fe;
M=fe.M;n=size(M,1);
assert(size(M,2)==n&&mod(n,6)==0,'BridgeWest:MassSchema','Expected six DOFs per node in the native full mass matrix.');
mass=zeros(3,1);
for j=1:3
    r=zeros(n,1);r(j:6:end)=1;mass(j)=full(r.'*M*r);
end
assert(all(isfinite(mass))&&all(mass>0),'BridgeWest:Mass','Invalid rigid-translation mass.');
assert((max(mass)-min(mass))/mean(mass)<1e-6,'BridgeWest:Mass','Directional translational masses disagree.');
components=table(["steel";"walkway_surface";"residual_added_mass"], ...
    [double(fe.steel_mass);double(fe.walkway_mass);double(fe.residual_mass)], ...
    ["depends on provisional geometry";"assigned surface density and width";"closes imposed total; not independently measured ballast"], ...
    'VariableNames',{'component','mass_kg','evidence_status'});
assert(abs(sum(components.mass_kg)-mean(mass))<1e-6*mean(mass),'BridgeWest:Mass','Mass partition sum does not match assembled mass.');
assert(isfield(fe.c,'mass_per_length_total_kg_m')&&isfield(fe.c,'length_m'),'BridgeWest:MassConfig','Missing nominal mass settings.');
imposed=double(fe.c.mass_per_length_total_kg_m)*double(fe.c.length_m);
rows=struct('event',{},'campaign',{},'inferred_total_kN',{},'weighed_total_kN',{}, ...
    'inferred_over_weighed',{},'weighed_over_inferred',{},'n_bogies_le182_kN',{});
events=bw_select_saved_events(outputRoot,'Event','all');allLoads=[];
for i=1:numel(events)
    eid=events{i};P=load(fullfile(outputRoot,'predictions',[eid '.mat']),'load_mean_kN');
    p=double(P.load_mean_kN(:));assert(numel(p)==8&&all(isfinite(p))&&all(p>0),'BridgeWest:Loads','Invalid eight-bogie load vector in %s.',eid);
    year=str2double(eid(1:4));W=bw_weighed_bogie_loads(year,o.WeightDir);tot=sum(p);allLoads=[allLoads;p]; %#ok<AGROW>
    rows(end+1)=struct('event',eid,'campaign',year,'inferred_total_kN',tot, ...
        'weighed_total_kN',W.total_load_kN,'inferred_over_weighed',tot/W.total_load_kN, ...
        'weighed_over_inferred',W.total_load_kN/tot,'n_bogies_le182_kN',sum(p<=182)); %#ok<AGROW>
end
T=struct2table(rows);
A=struct('mass_x_kg',mass(1),'mass_y_kg',mass(2),'mass_z_kg',mass(3), ...
    'configured_total_mass_kg',imposed,'mass_components',components,'train_load_audit',T, ...
    'n_passages',height(T),'n_bogie_means',numel(allLoads),'n_bogie_means_le182_kN',sum(allLoads<=182), ...
    'median_bogie_mean_kN',median(allLoads),'independent_bridge_mass_established',false, ...
    'bridge_mass_scale_applied',1,'prediction_files_modified',false);
writetable(components,fullfile(o.OutputDir,'bridge_mass_components.csv'));
writetable(T,fullfile(o.OutputDir,'train_load_audit.csv'));
save(fullfile(o.OutputDir,'mass_and_load_audit.mat'),'A');
fprintf('Native rigid-translation mass: %.3f tonnes. Configured total: %.3f tonnes.\n',mean(mass)/1000,imposed/1000);
fprintf('Residual assigned added mass: %.3f tonnes; not an independent ballast measurement.\n',fe.residual_mass/1000);
fprintf('%d of %d bogie posterior means are <=182 kN.\n',sum(allLoads<=182),numel(allLoads));
fprintf('No bridge-mass correction applied; no prediction files changed.\n');
disp(components);
end
