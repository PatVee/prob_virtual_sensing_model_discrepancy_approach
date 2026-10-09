function out=run_bridgewest_joint_physical_v1(packageRoot,v3ResultRoot,outputRoot,events)
%RUN_BRIDGEWEST_JOINT_PHYSICAL_V1 Joint model with shaker-updated physical M.
% The FE is physically updated BEFORE moving-load inference. Two shaker-only
% mass-matrix parameters are identified: dynamic participation of the residual
% track/ballast mass and its transverse redistribution. K is unchanged, M is
% updated, and the eigensystem is re-solved. Modal frequencies are therefore
% outputs of K*phi=lambda*M*phi; they are not overwritten for the baseline FE.
%
% Development / calibration information:
%   - shaker low-force f1/f2 + four-station shapes -> physical M update;
%   - 2022 independent static bogie weighing + A031/B053 -> one frozen
%     strain-recovery observation scale, retained from the joint model.
%
% 2025 A031/B053 infer passage loads. 2025 static weighing, target strains and
% PLS displacement are opened only after every requested prediction is frozen.
% No shaker modal-residue eigenvector scaling is applied in this branch.
scriptRoot=fileparts(mfilename('fullpath'));
if nargin<1||isempty(packageRoot),packageRoot=scriptRoot;end
packageRoot=char(packageRoot);
if nargin<2||isempty(v3ResultRoot)
    cand=fullfile(fileparts(packageRoot),'BridgeWest_Stabil_Shaker_v3_results');
    if isfolder(cand),v3ResultRoot=cand;else,error('Pass BridgeWest_Stabil_Shaker_v3_results as v3ResultRoot.');end
end
v3ResultRoot=char(v3ResultRoot);
if nargin<3||isempty(outputRoot),outputRoot=fullfile(scriptRoot,'joint_physical_results');end
if ~isfolder(outputRoot),mkdir(outputRoot);end

addpath(fullfile(packageRoot,'matlab'));
addpath(fullfile(packageRoot,'dependencies','Stabil'));
addpath(fullfile(packageRoot,'shaker_informed','matlab'));
addpath(fullfile(packageRoot,'production','matlab'));
addpath(fullfile(packageRoot,'production','matlab','raw_functions'));
sc=jsondecode(fileread(fullfile(packageRoot,'shaker_informed','discrepancy_config.json')));
c=jsondecode(fileread(fullfile(packageRoot,'production','config','prediction_config.json')));
evset=jsondecode(fileread(fullfile(packageRoot,'production','config','event_sets.json')));
copyFields={'cov_ridge','reference_noise_sigma_microstrain','target_noise_sigma_microstrain','displacement_noise_sigma_mm','interval_z'};
for k=1:numel(copyFields),sc.(copyFields{k})=c.(copyFields{k});end
if nargin<4||isempty(events),events=arrayfun(@(x)sprintf('2025_run%02d',x),double(evset.runs_2025(:)),'UniformOutput',false);else,events=cellstr(string(events));end

%% FE + shaker-only true physical mass-matrix update
feDir=fullfile(outputRoot,'fe_model');if ~isfolder(feDir),mkdir(feDir);end
feRaw=run_bridgewest_lorenzen(fullfile(packageRoot,'config','model.json'),feDir);
save(fullfile(feDir,'bridgewest_lorenzen_fe_raw.mat'),'feRaw','-v7.3');
[fe,physicalUpdate]=bw_shaker_physical_mass_update(feRaw,packageRoot,'Calibrate',true,'Verbose',true);
feDyn=fe;
sc.bending_mode_fe_index=physicalUpdate.bending_mode_fe_index;
sc.torsion_mode_fe_index=physicalUpdate.torsion_mode_fe_index;
save(fullfile(feDir,'bridgewest_lorenzen_fe.mat'),'fe','-v7.3');
modalTable=table([physicalUpdate.bending_mode_fe_index;physicalUpdate.torsion_mode_fe_index], ...
    physicalUpdate.updated_frequency_Hz(:),physicalUpdate.target_shaker_frequency_Hz(:), ...
    physicalUpdate.shaker_spatial_MAC(:),physicalUpdate.modal_residue_norm(:), ...
    physicalUpdate.shaker_modal_residue_norm(:),physicalUpdate.modal_residue_ratio_to_shaker(:), ...
    'VariableNames',{'FE_mode_index','physical_updated_frequency_Hz','shaker_frequency_Hz', ...
    'shaker_spatial_MAC','FE_residue_norm','shaker_residue_norm','residue_ratio_to_shaker'});
modalCompare=table({'bending';'torsion'}, ...
    [physicalUpdate.raw_bending_mode_fe_index;physicalUpdate.raw_torsion_mode_fe_index], ...
    physicalUpdate.raw_frequency_Hz(:),physicalUpdate.raw_shaker_spatial_MAC(:), ...
    physicalUpdate.raw_modal_residue_ratio_to_shaker(:), ...
    [physicalUpdate.bending_mode_fe_index;physicalUpdate.torsion_mode_fe_index], ...
    physicalUpdate.updated_frequency_Hz(:),physicalUpdate.shaker_spatial_MAC(:), ...
    physicalUpdate.modal_residue_ratio_to_shaker(:),physicalUpdate.target_shaker_frequency_Hz(:), ...
    'VariableNames',{'shaker_family','raw_FE_mode_index','raw_FE_frequency_Hz','raw_spatial_MAC', ...
    'raw_residue_ratio_to_shaker','updated_FE_mode_index','updated_FE_frequency_Hz', ...
    'updated_spatial_MAC','updated_residue_ratio_to_shaker','shaker_frequency_Hz'});
writetable(modalTable,fullfile(outputRoot,'physical_modal_calibration.csv'));
writetable(modalCompare,fullfile(outputRoot,'physical_modal_update_comparison.csv'));
writetable(modalTable,fullfile(outputRoot,'modal_residue_calibration.csv')); % plotter compatibility
writeJson(fullfile(outputRoot,'physical_update.json'),physicalUpdate);

%% 2022 load-informed strain-recovery calibration (A031/B053 only)
W=readtable(fullfile(packageRoot,'audit','verified_static_bogie_loads.csv'));
q22=sortrows(W(W.campaign==2022,:),'bogie').static_load_kN;assert(numel(q22)==8);q22N=double(q22(:))*1000;
refs=double(c.reference_model_indices(:)).';calRows=struct([]);
for ir=1:numel(evset.runs_2022)
    eid=sprintf('2022_run%02d',double(evset.runs_2022(ir)));pin=loadPermitted(v3ResultRoot,eid);
    B=bw_joint_nominal_basis(feDyn,pin.speed_m_s,pin.time_s,sc);
    pred=squeeze(sum(B{1}(:,refs,:).*reshape(q22N,1,1,8),3));y=pin.reference_strain_microstrain;
    sA=scalarFit(pred(:,1),y(:,1));sB=scalarFit(pred(:,2),y(:,2));sE=scalarFit(pred(:),y(:));
    calRows=appendStruct(calRows,struct('event',eid,'scale_event',sE,'scale_A031',sA,'scale_B053',sB));
end
calTable=struct2table(calRows);strainScale=median(calTable.scale_event,'omitnan');
writetable(calTable,fullfile(outputRoot,'2022_reference_scale_by_event.csv'));
fid=fopen(fullfile(outputRoot,'frozen_calibration.txt'),'w');fprintf(fid,'Physical residual participation %.15g\nPhysical outer redistribution %.15g\nUpdated bending mode %d, %.15g Hz\nUpdated torsion mode %d, %.15g Hz\n2022 median strain-recovery scale %.15g\n',physicalUpdate.residual_dynamic_participation,physicalUpdate.outer_redistribution_fraction,physicalUpdate.bending_mode_fe_index,physicalUpdate.updated_frequency_Hz(1),physicalUpdate.torsion_mode_fe_index,physicalUpdate.updated_frequency_Hz(2),strainScale);fclose(fid);

%% 2025 reference-only predictions; no 2025 validation source opened here
predDir=fullfile(outputRoot,'predictions');if ~isfolder(predDir),mkdir(predDir);end
loadRows=struct([]);
for ie=1:numel(events)
    eid=char(events{ie});fprintf('[%d/%d] %s\n',ie,numel(events),eid);pin=loadPermitted(v3ResultRoot,eid);
    [B0,components,meta]=bw_joint_basis_components(feDyn,pin.speed_m_s,pin.time_s,sc,strainScale);
    [pred,dg]=bw_predict_continuous_joint(B0,components,pin.reference_strain_microstrain,refs,c);
    q=double(pred.load_mean_kN(:));S=double(dg.load_posterior.samples);tot=sum(S,2);ql=empQuantile(tot,[.025 .5 .975]);
    provenance=struct('event',eid,'physical_update',physicalUpdate,'strain_recovery_scale',strainScale, ...
        'joint_meta',meta,'protected_targets_used',false,'PLS_used',false,'static_2025_used',false);
    save(fullfile(predDir,[eid '_prediction.mat']),'pred','provenance','-v7.3');
    row=struct('event',eid,'total_kN',sum(q),'total_q025_kN',ql(1),'total_q50_kN',ql(2),'total_q975_kN',ql(3), ...
        'n_le182',sum(q<=182),'n_ge358',sum(q>=358));
    for j=1:8,row.(sprintf('q%d_kN',j))=q(j);end
    loadRows=appendStruct(loadRows,row);
end
L=struct2table(loadRows);writetable(L,fullfile(outputRoot,'reference_only_2025_loads_UNSCORED.csv'));

%% Validation opened only after every prediction is frozen
q25=sortrows(W(W.campaign==2025,:),'bogie').static_load_kN;staticTotal=sum(q25);
L.static_total_kN=repmat(staticTotal,height(L),1);L.total_error_percent=100*(L.total_kN/staticTotal-1);
L.static_inside_95=(L.total_q025_kN<=staticTotal & L.total_q975_kN>=staticTotal);
writetable(L,fullfile(outputRoot,'2025_load_validation.csv'));
strainRows=struct([]);dispRows=struct([]);
for i=1:height(L)
    eid=char(string(L.event(i)));P=load(fullfile(predDir,[eid '_prediction.mat']),'pred');pred=P.pred;
    vf=fullfile(v3ResultRoot,'validation_targets',[eid '.mat']);if isfile(vf),V=load(vf);strainRows=appendStructArray(strainRows,scoreStrain(eid,pred,V,c));end
    pf=fullfile(v3ResultRoot,'displacement_validation',[eid '_measured_PLS_snapshot.csv']);if isfile(pf),T=readtable(pf,'TextType','string');dispRows=appendStructArray(dispRows,scoreDisp(eid,pred,T,feDyn));end
end
if isempty(strainRows),ST=table();else,ST=struct2table(strainRows);writetable(ST,fullfile(outputRoot,'2025_strain_validation.csv'));end
if isempty(dispRows),DT=table();else,DT=struct2table(dispRows);writetable(DT,fullfile(outputRoot,'2025_displacement_validation.csv'));end

summary=struct('physical_update',physicalUpdate,'strain_recovery_scale',strainScale, ...
    'median_abs_load_error_percent',median(abs(L.total_error_percent),'omitnan'), ...
    'static_inside_95_fraction',mean(L.static_inside_95));
if ~isempty(DT)
    for g={'BA','BB'}
        q=strcmpi(string(DT.girder),g{1});summary.([g{1} '_median_peak_ratio'])=median(DT.peak_ratio(q),'omitnan');summary.([g{1} '_median_nrmse'])=median(DT.spatial_nrmse(q),'omitnan');
    end
end
writeJson(fullfile(outputRoot,'summary.json'),summary);
out=struct('outputRoot',outputRoot,'modal_residue_calibration',modalTable,'modal_update_comparison',modalCompare,'strain_scale_calibration',calTable,'loads',L,'strain',ST,'displacement',DT,'summary',summary);
end

function pin=loadPermitted(root,eid)
f=fullfile(root,'permitted_inputs',[eid '.mat']);assert(isfile(f),'Missing %s',f);pin=load(f);pin.time_s=double(pin.time_s(:));pin.speed_m_s=double(pin.speed_m_s);pin.reference_strain_microstrain=double(pin.reference_strain_microstrain);
end
function s=scalarFit(x,y),x=double(x(:));y=double(y(:));ok=isfinite(x)&isfinite(y);x=x(ok);y=y(ok);s=(x.'*y)/max(x.'*x,eps);end
function v=modeAtTop(fe,g,x,j),xn=double(fe.x(:));k=find(xn<=x,1,'last');if isempty(k),k=1;end;k=min(k,numel(xn)-1);a=(x-xn(k))/max(xn(k+1)-xn(k),eps);ids=double(fe.top(g,[k k+1]));dof=6*(ids-1)+3;vv=double(fe.P(dof,j));v=(1-a)*vv(1)+a*vv(2);end
function r=residueNorm(fe,j),left=double(fe.c.bearing_x_m(1));xs=[left+4.1 left+8.2];resp=[modeAtTop(fe,1,xs(1),j);modeAtTop(fe,2,xs(1),j);modeAtTop(fe,1,xs(2),j);modeAtTop(fe,2,xs(2),j)];force=modeAtTop(fe,2,8.7,j);r=abs(force)*norm(resp);end
function q=empQuantile(x,p),x=sort(double(x(:)));n=numel(x);q=zeros(size(p));for i=1:numel(p),u=1+(n-1)*p(i);a=floor(u);b=ceil(u);if a==b,q(i)=x(a);else,q(i)=x(a)+(u-a)*(x(b)-x(a));end;end;end
function rows=scoreStrain(eid,p,V,c)
rows=struct([]);y=double(V.Y);flags=string(V.flags(:));for j=1:numel(c.target_model_indices),if flags(j)~="ok"||any(~isfinite(y(:,j))),continue;end;midx=double(c.target_model_indices(j));mu=double(p.strain_mean(:,midx));lo=double(p.strain_lo(:,midx));hi=double(p.strain_hi(:,midx));yy=y(:,j);rg=max(yy)-min(yy);C=corrcoef(mu,yy);if numel(C)>=4,rr=C(1,2);else,rr=NaN;end;rows=appendStruct(rows,struct('event',eid,'target',char(string(c.target_labels(j))),'nrmse',sqrt(mean((mu-yy).^2))/max(rg,eps),'correlation',rr,'range_ratio',(max(mu)-min(mu))/max(rg,eps),'coverage',mean(yy>=lo&yy<=hi)));end
end
function rows=scoreDisp(eid,p,T,fe)
rows=struct([]);x=double(fe.observation_maps.x_dataset_m(:));for gg={'BA','BB'},g=gg{1};if strcmp(g,'BA'),Y=double(p.disp_BA_mean);SD=double(p.disp_BA_sd);else,Y=double(p.disp_BB_mean);SD=double(p.disp_BB_sd);end;[xs,ord]=sort(x);mid=mean([min(x) max(x)]);[~,jm]=min(abs(x-mid));[~,it]=max(abs(Y(:,jm)));mu=Y(it,ord).';sd=SD(it,ord).';q=T(strcmpi(string(T.axis),g),:);if height(q)<3,continue;end;if ismember('x_m',q.Properties.VariableNames),xm=double(q.x_m);ym=double(q.measured_mm);else,xm=double(q.x);ym=double(q.measured);end;[xm,io]=sort(xm);ym=ym(io);yp=interp1(xs,mu,xm,'linear',NaN);sp=interp1(xs,sd,xm,'linear',NaN);ok=isfinite(ym)&isfinite(yp);ym=ym(ok);yp=yp(ok);sp=sp(ok);rg=max(ym)-min(ym);C=corrcoef(yp,ym);if numel(C)>=4,rr=C(1,2);else,rr=NaN;end;lo=yp-1.959963984540054*sp;hi=yp+1.959963984540054*sp;rows=appendStruct(rows,struct('event',eid,'girder',g,'peak_ratio',max(abs(yp))/max(max(abs(ym)),eps),'spatial_nrmse',sqrt(mean((yp-ym).^2))/max(rg,eps),'correlation',rr,'profile_coverage',mean(ym>=lo&ym<=hi)));end
end
function A=appendStruct(A,row),if isempty(A),A=row;else,A(end+1)=row;end,end
function A=appendStructArray(A,B),if isempty(B),return;end;if isempty(A),A=B;else,A=[A B];end,end
function writeJson(file,x),fid=fopen(file,'w');assert(fid>=0);c=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(x,'PrettyPrint',true));end
