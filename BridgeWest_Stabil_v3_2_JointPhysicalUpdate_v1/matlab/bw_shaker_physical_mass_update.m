function [fe,upd]=bw_shaker_physical_mass_update(feRaw,packageRoot,varargin)
%BW_SHAKER_PHYSICAL_MASS_UPDATE Re-identify M and re-solve the FE eigensystem.
%
% This is a true matrix update: the physical FE mass matrix is modified first,
% then K*phi=lambda*M*phi is solved again. Frequencies are NOT overwritten.
%
% Two shaker-only parameters are identified:
%   eta    : dynamic participation fraction of the residual nonstructural
%            track/ballast mass already present in feRaw.added_mass;
%   lambda : fraction of participating residual mass redistributed, at each
%            longitudinal x-slice, from the tributary deck distribution to
%            the two outer girder/deck lines. Total mass in each x-slice is
%            preserved before eta is applied.
%
% The sourced static mass inventory is not changed. eta is an equivalent
% dynamic participation factor for the shaker-frequency model.
%
% No train-passage strain, static weighing, target strain or PLS data are read.

ip=inputParser;
addParameter(ip,'Calibrate',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Verbose',true,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});opt=ip.Results;

cfg=jsondecode(fileread(fullfile(packageRoot,'config','physical_update.json')));
ens=bw_shaker_ensemble(fullfile(packageRoot,'shaker_informed'));
shapeFile=fullfile(packageRoot,'shaker_informed','data','stw_spatial_mode_shapes_2025.csv');
S=readtable(shapeFile,'TextType','string');
obs1=double(S.mode1_STW(:));obs2=double(S.mode2_STW(:));
target=[double(ens.baseline.f1_Hz);double(ens.baseline.f2_Hz)];

% Split the already assembled added mass into walkway and residual mass.
[nodalWalk,resUniform,resOuter]=massParts(feRaw);
nominalTotal=double(feRaw.c.mass_per_length_total_kg_m*feRaw.c.length_m);
staticResidual=double(feRaw.residual_mass);
staticSteel=double(feRaw.steel_mass);staticWalk=double(feRaw.walkway_mass);

etaB=double(cfg.residual_dynamic_participation_bounds(:));
lamB=double(cfg.outer_redistribution_bounds(:));
eta0=double(cfg.initial_residual_dynamic_participation);
lam0=double(cfg.initial_outer_redistribution_fraction);
z0=[invBound(eta0,etaB);invBound(lam0,lamB)];

rawNM=feRaw.c.number_of_modes;
calNM=min(rawNM,double(cfg.calibration_number_of_modes));
objTrace=[];
if opt.Calibrate
    oo=optimset('Display','off','MaxFunEvals',double(cfg.max_function_evaluations), ...
        'MaxIter',double(cfg.max_function_evaluations),'TolX',2e-5,'TolFun',1e-10);
    z=fminsearch(@objective,z0,oo);
else
    z=z0;
end
eta=bound(z(1),etaB);lambda=bound(z(2),lamB);

% Raw and final full-order modal mappings.  The raw mapping is retained for
% an explicit pre-/post-update audit; the final mapping is re-identified after
% M is changed because modal order is allowed to change.
rawMap=mapModes(feRaw,S,obs1,obs2,cfg);
rawIb=rawMap.bending_index;rawIt=rawMap.torsion_index;
rawRes=[residueNorm(feRaw,rawIb);residueNorm(feRaw,rawIt)];

fe=applyUpdate(feRaw,eta,lambda,nodalWalk,resUniform,resOuter,rawNM);
map=mapModes(fe,S,obs1,obs2,cfg);
ib=map.bending_index;it=map.torsion_index;
res=[residueNorm(fe,ib);residueNorm(fe,it)];
OLD_RES1=2.142090694997057e-05;OLD_RES2=5.204858470581088e-05;
OLD_TO_MEAS1=1.59272552937298;OLD_TO_MEAS2=0.945167979654921;
measRes=[OLD_RES1*OLD_TO_MEAS1;OLD_RES2*OLD_TO_MEAS2];

effectiveResidual=eta*staticResidual;
effectiveTotal=staticSteel+staticWalk+effectiveResidual;
upd=struct();
upd.residual_dynamic_participation=eta;
upd.outer_redistribution_fraction=lambda;
upd.nominal_static_total_mass_kg=nominalTotal;
upd.static_residual_mass_kg=staticResidual;
upd.effective_dynamic_residual_mass_kg=effectiveResidual;
upd.effective_dynamic_total_mass_kg=effectiveTotal;
upd.effective_dynamic_mass_fraction_of_nominal=effectiveTotal/nominalTotal;
upd.raw_bending_mode_fe_index=rawIb;
upd.raw_torsion_mode_fe_index=rawIt;
upd.raw_frequency_Hz=[double(feRaw.freq(rawIb));double(feRaw.freq(rawIt))];
upd.raw_shaker_spatial_MAC=[rawMap.bending_MAC;rawMap.torsion_MAC];
upd.raw_top_line_participation=[rawMap.bending_top_ratio;rawMap.torsion_top_ratio];
upd.raw_modal_residue_norm=rawRes;
upd.raw_modal_residue_ratio_to_shaker=rawRes./measRes;
upd.bending_mode_fe_index=ib;
upd.torsion_mode_fe_index=it;
upd.updated_frequency_Hz=[double(fe.freq(ib));double(fe.freq(it))];
upd.target_shaker_frequency_Hz=target;
upd.frequency_error_percent=100*(upd.updated_frequency_Hz./target-1);
upd.shaker_spatial_MAC=[map.bending_MAC;map.torsion_MAC];
upd.top_line_participation=[map.bending_top_ratio;map.torsion_top_ratio];
upd.modal_residue_norm=res;
upd.shaker_modal_residue_norm=measRes;
upd.modal_residue_ratio_to_shaker=res./measRes;
upd.calibration_uses_train_passages=false;
upd.objective_trace=objTrace;
fe.physical_update=upd;
fe.observation_maps=bw_lorenzen_observation_maps(fe);

if opt.Verbose
    fprintf('\nShaker-informed physical mass update\n');
    fprintf('  residual dynamic participation eta = %.8f\n',eta);
    fprintf('  outer redistribution lambda        = %.8f\n',lambda);
    fprintf('  effective dynamic total mass       = %.3f kg (%.3f%% nominal)\n',effectiveTotal,100*effectiveTotal/nominalTotal);
    fprintf('  bending mode %d: %.6f Hz; MAC %.6f; residue ratio %.4f\n',ib,fe.freq(ib),map.bending_MAC,res(1)/measRes(1));
    fprintf('  torsion mode %d: %.6f Hz; MAC %.6f; residue ratio %.4f\n',it,fe.freq(it),map.torsion_MAC,res(2)/measRes(2));
end

    function val=objective(zz)
        et=bound(zz(1),etaB);la=bound(zz(2),lamB);
        ftmp=applyUpdate(feRaw,et,la,nodalWalk,resUniform,resOuter,calNM);
        mm=mapModes(ftmp,S,obs1,obs2,cfg);
        ff=[double(ftmp.freq(mm.bending_index));double(ftmp.freq(mm.torsion_index))];
        rf=log(ff./target);
        rs=[1-mm.bending_MAC;1-mm.torsion_MAC];
        val=double(cfg.objective_frequency_log_weight)*sum(rf.^2)+ ...
            double(cfg.objective_shape_weight)*sum(rs.^2);
        objTrace(end+1,:)=[et la ff(:).' mm.bending_MAC mm.torsion_MAC val]; %#ok<AGROW>
    end
end

function fe=applyUpdate(fe0,eta,lambda,walk,resU,resO,nm)
fe=fe0;n=size(fe.M,1);old=double(fe.added_mass(:));new=walk+eta*((1-lambda)*resU+lambda*resO);
delta=new-old;dd=zeros(n,1);
for j=1:3,dd(j:6:end)=delta;end
fe.M=fe.M+spdiags(dd,0,n,n);fe.added_mass=new;
fe.effective_residual_mass=eta*double(fe.residual_mass);
fe.dynamic_mass_participation=eta;fe.outer_mass_fraction=lambda;
oldNM=fe.c.number_of_modes;fe.c.number_of_modes=nm;fe=bw_lorenzen_modes(fe);fe.c.number_of_modes=oldNM;
end

function [walk,resU,resO]=massParts(fe)
nn=size(fe.nodes,1);walkArea=zeros(nn,1);c=fe.c;
for e=1:size(fe.beams,1)
    row=fe.beams(e,:);ids=row(1:2);p=fe.nodes(ids,:);
    if row(7)==2 && abs(p(2,1)-p(1,1))>1e-7
        le=norm(p(2,:)-p(1,:));walkArea(ids)=walkArea(ids)+le*c.walkway_width_m/2;
    end
end
walk=walkArea*c.walkway_surface_mass_kg_m2;
resU=double(fe.added_mass(:))-walk;
assert(min(resU)>-1e-7,'Residual nodal mass became negative.');resU=max(resU,0);
resO=zeros(nn,1);xu=unique(round(fe.nodes(:,1),9));
for ix=1:numel(xu)
    ids=find(abs(fe.nodes(:,1)-xu(ix))<5e-9);mx=sum(resU(ids));if mx<=0,continue;end
    edge=zeros(2,1);
    for g=1:2
        target=[xu(ix) c.girder_y_m(g) c.deck_z_m];
        [~,k]=min(sum((fe.nodes(ids,:)-target).^2,2));edge(g)=ids(k);
    end
    resO(edge)=resO(edge)+mx/2;
end
assert(abs(sum(resO)-sum(resU))<1e-6*max(sum(resU),1),'Outer redistribution does not conserve residual mass.');
end

function m=mapModes(fe,S,obs1,obs2,cfg)
nm=numel(fe.freq);mac1=-inf(nm,1);mac2=-inf(nm,1);topRatio=zeros(nm,1);
topDof=6*(double(fe.top(:))-1)+3;allZDof=(3:6:size(fe.P,1))';
for j=1:nm
    if fe.freq(j)>double(cfg.mode_mapping_max_frequency_Hz),continue;end
    pred=zeros(height(S),1);
    for k=1:height(S)
        g=1+strcmpi(char(S.girder(k)),'B');
        pred(k)=modeAtTop(fe,g,double(fe.c.bearing_x_m(1))+double(S.x_FE_m(k)),j);
    end
    mac1(j)=MAC(pred,obs1);mac2(j)=MAC(pred,obs2);
    topRatio(j)=max(abs(fe.P(topDof,j)))/max(max(abs(fe.P(allZDof,j))),eps);
end
valid=topRatio>=double(cfg.minimum_top_line_participation) & isfinite(mac1);
assert(any(valid),'No physically visible FE modes available for shaker mapping.');
a=mac1;a(~valid)=-inf;[bMAC,ib]=max(a);
a=mac2;a(~valid)=-inf;a(ib)=-inf;[tMAC,it]=max(a);
assert(isfinite(bMAC)&&isfinite(tMAC),'Could not map two distinct shaker modes.');
m=struct('bending_index',ib,'torsion_index',it,'bending_MAC',bMAC,'torsion_MAC',tMAC, ...
    'bending_top_ratio',topRatio(ib),'torsion_top_ratio',topRatio(it));
end

function v=modeAtTop(fe,g,x,j)
xn=double(fe.x(:));k=find(xn<=x,1,'last');if isempty(k),k=1;end;k=min(k,numel(xn)-1);
a=(x-xn(k))/max(xn(k+1)-xn(k),eps);ids=double(fe.top(g,[k k+1]));dof=6*(ids-1)+3;vv=double(fe.P(dof,j));v=(1-a)*vv(1)+a*vv(2);
end
function r=residueNorm(fe,j)
left=double(fe.c.bearing_x_m(1));xs=[left+4.1 left+8.2];resp=[modeAtTop(fe,1,xs(1),j);modeAtTop(fe,2,xs(1),j);modeAtTop(fe,1,xs(2),j);modeAtTop(fe,2,xs(2),j)];force=modeAtTop(fe,2,8.7,j);r=abs(force)*norm(resp);
end
function v=MAC(a,b),a=double(a(:));b=double(b(:));v=abs(a.'*b)^2/max((a.'*a)*(b.'*b),eps);end
function x=bound(z,b),x=b(1)+(b(2)-b(1))/(1+exp(-z));end
function z=invBound(x,b),q=(x-b(1))/(b(2)-b(1));q=min(max(q,1e-8),1-1e-8);z=log(q/(1-q));end
