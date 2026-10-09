function H=kw51_basis(B,E,t,kin,outputs)
% Unit axle influence coefficients, multiplied by 100 kN.
if nargin<5,outputs=1:7;end
t=t(:);tau=E.rail_offsets_s(:).';v=kin(1);entry=kin(2);direction=kin(3);acc=kin(4);
dt=t-entry;vr=v+acc*(E.rail_first_time_s-entry);
x=v*dt+.5*acc*dt.^2-(vr*tau+.5*acc*tau.^2);
if direction<0,x=116.316-x;end
speed2=(v+acc*dt).^2;tr=1+strcmp(E.track,'B');
grid=B.xgrid(:);dx=grid(2)-grid(1);assert(max(abs(diff(grid)-dx))<1e-10,'FE grid must be uniform.');
valid=x>=grid(1)&x<=grid(end);q=(x-grid(1))/dx;
i=min(max(floor(q)+1,1),numel(grid)-1);a=q-(i-1);
H=zeros(numel(t),numel(outputs),numel(tau));
for j=1:numel(outputs)
 h0=B.H0(:,outputs(j),tr);h2=B.H2(:,outputs(j),tr);
 u0=reshape(h0(i(:)),size(x));d0=reshape(h0(i(:)+1)-h0(i(:)),size(x));
 u2=reshape(h2(i(:)),size(x));d2=reshape(h2(i(:)+1)-h2(i(:)),size(x));
 z=((u0+a.*d0)+speed2.*(u2+a.*d2)).*valid;
 H(:,j,:)=reshape(z*1e5,numel(t),1,numel(tau));
end
end
