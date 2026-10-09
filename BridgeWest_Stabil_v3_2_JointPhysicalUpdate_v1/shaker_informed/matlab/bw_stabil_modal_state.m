function state = bw_stabil_modal_state(fe,shaker,c)
%BW_STABIL_MODAL_STATE Apply one shaker state to the detailed Stabil modal ROM.
% Geometry, mass, mode shapes and the documented top-flange restraint remain
% fixed. Only the identified first global bending and torsional eigenvalues
% and their damping anchors are replaced by shaker observations.
if nargin<3 || isempty(c)
    root=fileparts(fileparts(mfilename('fullpath')));
    c=jsondecode(fileread(fullfile(root,'discrepancy_config.json')));
end
ib=c.bending_mode_fe_index; it=c.torsion_mode_fe_index;
freq=fe.freq(:);
assert(ib>=1 && ib<=numel(freq) && it>=1 && it<=numel(freq),'Configured modal indices exceed FE modes.');
freq(ib)=double(shaker.f1_Hz); freq(it)=double(shaker.f2_Hz);
omega=2*pi*freq;
% Rayleigh damping is anchored to the two shaker-identified global modes.
A=[1/(2*omega(ib)) omega(ib)/2; 1/(2*omega(it)) omega(it)/2];
ab=A\[double(shaker.zeta1);double(shaker.zeta2)];
zeta=ab(1)./(2*omega)+ab(2).*omega/2;
zeta=max(min(zeta,0.08),0.002);
zeta(ib)=double(shaker.zeta1); zeta(it)=double(shaker.zeta2);
sel=find(fe.freq(:)<=c.dynamic_max_frequency_Hz);
sel=unique([sel(:);ib;it]);
state=struct('freq_Hz',freq,'omega_rad_s',omega,'zeta',zeta,'rayleigh_ab',ab, ...
    'selected_modes',sel,'bending_index',ib,'torsion_index',it, ...
    'f1_Hz',double(shaker.f1_Hz),'f2_Hz',double(shaker.f2_Hz), ...
    'zeta1',double(shaker.zeta1),'zeta2',double(shaker.zeta2));
end
