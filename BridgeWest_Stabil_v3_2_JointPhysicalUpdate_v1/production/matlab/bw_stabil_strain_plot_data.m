function D=bw_stabil_strain_plot_data(outputRoot,eid,dataRoot)
%BW_STABIL_STRAIN_PLOT_DATA Read saved predictions and assessment-only targets.
% No fitting, rescaling, load inversion, or change to prediction files occurs.
if nargin<3,dataRoot='';end
eid=char(string(eid));pf=fullfile(outputRoot,'predictions',[eid '.mat']);
assert(isfile(pf),'BridgeWest:MissingPrediction','Prediction file is absent: %s',pf);
P=load(pf);assert(strcmp(char(string(P.event)),eid),'BridgeWest:EventMismatch','Prediction event ID mismatch.');
t=double(P.time_s(:));
assert(numel(t)>1&&all(isfinite(t))&&all(diff(t)>0),'BridgeWest:Time','Invalid prediction time vector.');
labels={'A058','A070','B031','B081'};
measured={'str_x_BA_uf_058','str_x_BA_uf_070','str_x_BB_uf_031','str_x_BB_uf_081'};
assert(isfield(P,'model_channel_order'),'BridgeWest:Schema','Prediction has no model_channel_order.');
order=cellstr(string(P.model_channel_order));idx=zeros(1,4);
for j=1:4
    q=find(strcmp(order,labels{j}));
    assert(numel(q)==1,'BridgeWest:Schema','Expected exactly one model channel %s in %s.',labels{j},eid);idx(j)=q;
end
Y=double(P.strain_mean_microstrain);L=double(P.strain_lo_microstrain);U=double(P.strain_hi_microstrain);
assert(isequal(size(Y),size(L),size(U))&&size(Y,1)==numel(t)&&size(Y,2)>=max(idx), ...
    'BridgeWest:Schema','Strain array dimensions are inconsistent in %s.',eid);
Y=Y(:,idx);L=L(:,idx);U=U(:,idx);
assert(all(isfinite(Y(:)))&&all(isfinite(L(:)))&&all(isfinite(U(:)))&&all(L(:)<=U(:)), ...
    'BridgeWest:Schema','Non-finite predictions or invalid intervals in %s.',eid);
D=struct('event',eid,'time_s',t,'labels',{labels},'measured_channels',{measured}, ...
    'mean',Y,'lo',L,'hi',U,'measurement',nan(numel(t),4), ...
    'flags',{{'unavailable','unavailable','unavailable','unavailable'}}, ...
    'available',false(1,4),'primary_score',false(1,4),'nrmse',nan(1,4), ...
    'correlation',nan(1,4),'coverage',nan(1,4),'measurement_source','not available');
zf=fullfile(outputRoot,'validation_targets',[eid '.mat']);
if isfile(zf)
    Z=load(zf);D.measurement_source=zf;
elseif ~isempty(dataRoot)
    % Only this selected passage is read, after its prediction already exists.
    here=fileparts(mfilename('fullpath'));addpath(here,fullfile(here,'raw_functions'));
    sf=fullfile(outputRoot,'run_settings.json');assert(isfile(sf),'BridgeWest:Settings','Missing run_settings.json.');
    s=jsondecode(fileread(sf));
    assert(isfield(s,'preprocessing_config'),'BridgeWest:Settings','Saved settings lack preprocessing_config.');
    A=bw_prediction_zip_index(char(dataRoot),fullfile(outputRoot,'validation_archive_cache'));
    Z=bw_prepare_validation_targets(A,eid,s.preprocessing_config,outputRoot,measured,P.speed_m_s,P.entry_s,t);
    D.measurement_source=zf;
else
    return;
end
assert(isfield(Z,'channels')&&isfield(Z,'flags')&&isfield(Z,'Y'),'BridgeWest:ValidationSchema','Invalid validation file for %s.',eid);
if isfield(Z,'speed')&&isfield(P,'speed_m_s')
    assert(abs(double(Z.speed)-double(P.speed_m_s))<1e-8,'BridgeWest:ValidationKinematics','Cached speed differs from frozen prediction for %s.',eid);
end
if isfield(Z,'entry')&&isfield(P,'entry_s')
    assert(abs(double(Z.entry)-double(P.entry_s))<1e-8,'BridgeWest:ValidationKinematics','Cached entry time differs from frozen prediction for %s.',eid);
end
zt=double(Z.time(:));tol=1e-8*max(1,max(abs(t)));
assert(numel(zt)==numel(t)&&all(abs(zt-t)<=tol),'BridgeWest:TimeMismatch', ...
    'Validation/prediction time grids differ in %s. Regenerate validation with the frozen prediction times.',eid);
chs=cellstr(string(Z.channels));fl=cellstr(string(Z.flags));
assert(size(Z.Y,1)==numel(t)&&size(Z.Y,2)==numel(chs)&&numel(fl)==numel(chs), ...
    'BridgeWest:ValidationSchema','Validation arrays/labels disagree for %s.',eid);
for j=1:4
    col=find(strcmp(chs,measured{j}));
    if isempty(col),continue;end
    assert(numel(col)==1,'BridgeWest:ValidationSchema','Duplicate validation channel %s.',measured{j});
    y=double(Z.Y(:,col));D.measurement(:,j)=y;D.flags{j}=fl{col};
    D.available(j)=all(isfinite(y));D.primary_score(j)=D.available(j)&&strcmpi(fl{col},'ok');
    if ~D.primary_score(j),continue;end
    rg=max(y)-min(y);
    if rg>0,D.nrmse(j)=sqrt(mean((Y(:,j)-y).^2))/rg;end
    if std(y)>0&&std(Y(:,j))>0,q=corrcoef(y,Y(:,j));D.correlation(j)=q(1,2);end
    D.coverage(j)=mean(y>=L(:,j)&y<=U(:,j));
end
end
