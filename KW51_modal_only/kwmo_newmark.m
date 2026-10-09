function q=kwmo_newmark(force,dt,frequency_Hz,zeta)
% Independent unit-modal-mass oscillators; average-acceleration Newmark.
% Zero displacement/velocity at first padded time sample.
w=2*pi*frequency_Hz(:).';zeta=zeta(:).';
assert(size(force,2)==numel(w)&&numel(zeta)==numel(w));
k=w.^2;c=2*zeta.*w;a0=4/dt^2;a1=2/dt;
ke=k+a0+a1*c;
n=size(force,1);q=zeros(size(force));v=zeros(1,numel(w));a=force(1,:);
for i=2:n
 old=q(i-1,:);
 rhs=force(i,:)+a0*old+4/dt*v+a+c.*(a1*old+v);
 q(i,:)=rhs./ke;
 anew=a0*(q(i,:)-old)-4/dt*v-a;
 v=v+dt/2*(a+anew);a=anew;
end
end
