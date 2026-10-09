function out=bw_stabil_paired_discrepancy(fe,v,tobs,loadsN,c,direction,saveFile)
%BW_STABIL_PAIRED_DISCREPANCY Build shaker-informed unit-bogie discrepancies.
% The primary result is out.delta_basis: M x T x O x 8 unit-bogie response
% differences. This is structural/modal discrepancy before assigning uncertain
% train loads. If loadsN is supplied, out.delta is also returned after load
% projection. A single 1x8 load vector is applied to every state; Mx8 permits
% an explicit quadrature design. No random load draw is generated internally.
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<5 || isempty(c),c=jsondecode(fileread(fullfile(root,'discrepancy_config.json')));end
if nargin<6 || isempty(direction),direction=1;end
if nargin<7,saveFile='';end
ens=bw_shaker_ensemble(root);S=ens.states;M=height(S);
base=ens.baseline;[B0,meta0]=bw_stabil_basis_shaker(fe,v,base,tobs,c,direction);
deltaBasis={zeros(M,numel(tobs),7,8),zeros(M,numel(tobs),41,8),zeros(M,numel(tobs),41,8)};
for m=1:M
    sh=struct('f1_Hz',S.f1_Hz(m),'f2_Hz',S.f2_Hz(m),'zeta1',S.zeta1(m),'zeta2',S.zeta2(m));
    Bm=bw_stabil_basis_shaker(fe,v,sh,tobs,c,direction);
    for k=1:3
        deltaBasis{k}(m,:,:,:)=reshape(Bm{k}-B0{k},1,size(B0{k},1),size(B0{k},2),8);
    end
end
out=struct('delta_basis',{deltaBasis},'state_table',S,'baseline',base,'nominal_basis',{B0}, ...
    'weights',double(S.weight),'development_years',[2024 2025],'nominal_meta',meta0, ...
    'note','Shaker discrepancy is stored before load projection. 2024 uses genuine same-run pairs; 2025 uses within-level independent marginal products.');
if nargin>=4 && ~isempty(loadsN)
    out.loads_N=double(loadsN);out.delta=bw_project_discrepancy_basis(deltaBasis,loadsN);
else
    out.loads_N=[];out.delta=[];
end
if ~isempty(saveFile),save(saveFile,'out','-v7.3');end
end
