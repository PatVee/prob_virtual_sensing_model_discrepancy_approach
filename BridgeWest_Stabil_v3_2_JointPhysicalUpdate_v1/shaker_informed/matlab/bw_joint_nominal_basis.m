function [B0,meta]=bw_joint_nominal_basis(fe,v,tobs,sc)
%BW_JOINT_NOMINAL_BASIS Weighted mean of three distributed two-rail bases.
root=fileparts(fileparts(mfilename('fullpath')));ens=bw_shaker_ensemble(root);
hx=[0.6 1.2 2.4];hy=[0.30 0.45 0.60];w=ones(3,1)/3;Bstate=cell(3,1);
for i=1:3
    c=sc;c.load_transfer_model='two_rail_distributed';c.rail_gauge_m=1.435;
    c.rail_spread_longitudinal_half_length_m=hx(i);c.rail_spread_transverse_half_width_m=hy(i);c.conserve_wheel_resultant=true;
    Bstate{i}=bw_stabil_basis_shaker(fe,v,ens.baseline,tobs,c,1);
end
B0=cell(1,3);
for k=1:3
    B0{k}=zeros(size(Bstate{1}{k}));
    for i=1:3,B0{k}=B0{k}+w(i)*Bstate{i}{k};end
end
meta=struct('hx_m',hx,'hy_m',hy,'weights',w,'baseline',ens.baseline);
end
