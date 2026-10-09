function [B0,components,meta]=bw_joint_basis_components(fe,v,tobs,sc,strainScale)
%BW_JOINT_BASIS_COMPONENTS Construct four independent Bridge West layers.
%   1 shaker frequency/damping states
%   2 two-rail load-transfer footprint states
%   3 common compliance response-gain states
%   4 strain-recovery response-gain states
if nargin<5||isempty(strainScale),strainScale=1;end
root=fileparts(fileparts(mfilename('fullpath')));ens=bw_shaker_ensemble(root);
% Three distributed two-rail transfer states. These are an epistemic
% sensitivity mixture; target strain, PLS and 2025 weighing do not set them.
hx=[0.6 1.2 2.4];hy=[0.30 0.45 0.60];wt=ones(3,1)/3;
pairs=cell(3,1);
for il=1:3
    cl=sc;cl.load_transfer_model='two_rail_distributed';cl.rail_gauge_m=1.435;
    cl.rail_spread_longitudinal_half_length_m=hx(il);cl.rail_spread_transverse_half_width_m=hy(il);
    cl.conserve_wheel_resultant=true;
    pairs{il}=bw_stabil_paired_discrepancy(fe,v,tobs,[],cl,1,'');
end
B0=cell(1,3);dyn=cell(1,3);tr=cell(1,3);
M=size(pairs{1}.delta_basis{1},1);
for k=1:3
    B0{k}=zeros(size(pairs{1}.nominal_basis{k}));dyn{k}=zeros(size(pairs{1}.delta_basis{k}));
    for il=1:3
        B0{k}=B0{k}+wt(il)*pairs{il}.nominal_basis{k};
        dyn{k}=dyn{k}+wt(il)*pairs{il}.delta_basis{k};
    end
    tr{k}=zeros(3,size(B0{k},1),size(B0{k},2),8);
    for il=1:3
        tr{k}(il,:,:,:)=reshape(pairs{il}.nominal_basis{k}-B0{k},1,size(B0{k},1),size(B0{k},2),8);
    end
end
% Apply frozen strain-recovery scale to every strain operator and strain delta.
B0{1}=strainScale*B0{1};dyn{1}=strainScale*dyn{1};tr{1}=strainScale*tr{1};
% Residual common-compliance uncertainty: response gains 0.8/1.0/1.2,
% centred so it widens the model without fitting a validation mean.
cg=[0.8 1.0 1.2];cw=[0.25 0.50 0.25].';comp=cell(1,3);
for k=1:3
    comp{k}=zeros(3,size(B0{k},1),size(B0{k},2),8);
    for i=1:3,comp{k}(i,:,:,:)=reshape((cg(i)-1)*B0{k},1,size(B0{k},1),size(B0{k},2),8);end
end
% Residual strain-recovery uncertainty: +/-8%; displacement unaffected.
rg=[0.92 1.00 1.08];rw=[0.25 0.50 0.25].';rec=cell(1,3);
rec{1}=zeros(3,size(B0{1},1),size(B0{1},2),8);
for i=1:3,rec{1}(i,:,:,:)=reshape((rg(i)-1)*B0{1},1,size(B0{1},1),size(B0{1},2),8);end
rec{2}=zeros(3,size(B0{2},1),size(B0{2},2),8);rec{3}=zeros(3,size(B0{3},1),size(B0{3},2),8);
components=struct('name',{},'delta_basis',{},'weights',{});
components(1)=struct('name','shaker_dynamic','delta_basis',{dyn},'weights',double(ens.states.weight));
components(2)=struct('name','load_transfer','delta_basis',{tr},'weights',wt);
components(3)=struct('name','common_compliance','delta_basis',{comp},'weights',cw);
components(4)=struct('name','strain_recovery','delta_basis',{rec},'weights',rw);
meta=struct('transfer_hx_m',hx,'transfer_hy_m',hy,'transfer_weights',wt, ...
    'compliance_gains',cg,'compliance_weights',cw,'recovery_gains',rg,'recovery_weights',rw, ...
    'strain_recovery_scale',strainScale,'shaker_state_count',M);
end
