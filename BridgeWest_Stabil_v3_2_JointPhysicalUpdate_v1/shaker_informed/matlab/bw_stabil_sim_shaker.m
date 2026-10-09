function [t,H,state]=bw_stabil_sim_shaker(fe,v,shaker,c,direction)
%BW_STABIL_SIM_SHAKER Modal moving-load simulation for one shaker state.
% H is nTime x nSelectedModes x 8 unit-bogie bases.
if nargin<5 || isempty(direction),direction=1;end
state=bw_stabil_modal_state(fe,shaker,c);sel=state.selected_modes;
span=diff(fe.c.bearing_x_m);dt=c.simulation_dt_s;
t=(-c.pre_entry_s:dt:(span+max(c.axle_offsets_m))/v+c.post_exit_s).';
omega=state.omega_rad_s(sel);kk=omega.^2;cc=2*state.zeta(sel).*omega;
P=fe.P(:,sel);nm=numel(sel);nb=8;
beta=.25;gamma=.5;a0=1/(beta*dt^2);a1=gamma/(beta*dt);a2=1/(beta*dt);a3=1/(2*beta)-1;
a4=gamma/beta-1;a5=dt*(gamma/(2*beta)-1);den=kk+a0+a1*cc;
q=zeros(nm,nb);qd=q;qdd=q;H=zeros(numel(t),nm,nb);
for itime=1:numel(t)
    F=bw_stabil_force_map(fe,t(itime),v,c,direction);
    pm=P.'*F;
    if itime==1
        qdd=pm-cc.*qd-kk.*q;
    else
        pe=pm+a0*q+a2*qd+a3*qdd+cc.*(a1*q+a4*qd+a5*qdd);
        qn=pe./den;an=a0*(qn-q)-a2*qd-a3*qdd;
        vn=qd+dt*((1-gamma)*qdd+gamma*an);
        q=qn;qd=vn;qdd=an;
    end
    H(itime,:,:)=reshape(q,1,nm,nb);
end
end
