function [q,best,iterations]=kw51_nelder_mead(fun,q,lower,upper,maxiter)
% Deterministic bounded Nelder-Mead, matching the study's simplex/search rules.
n=numel(q);sim=repmat(q(:).',n+1,1);
for j=1:n,if q(j)~=0,sim(j+1,j)=1.05*q(j);else,sim(j+1,j)=.00025;end;end
over=sim>upper;reflected=2*upper-sim;sim(over)=reflected(over);
sim=min(max(sim,lower),upper);f=zeros(n+1,1);
for j=1:n+1,f(j)=fun(sim(j,:));end
[f,ix]=sort(f);sim=sim(ix,:);iterations=1;
while iterations<maxiter
 if max(max(abs(sim(2:end,:)-sim(1,:))))<=.005 && max(abs(f(2:end)-f(1)))<=1e-5,break;end
 center=mean(sim(1:n,:),1);xr=min(max(2*center-sim(end,:),lower),upper);fr=fun(xr);shrink=false;
 if fr<f(1)
  xe=min(max(3*center-2*sim(end,:),lower),upper);fe=fun(xe);
  if fe<fr,sim(end,:)=xe;f(end)=fe;else,sim(end,:)=xr;f(end)=fr;end
 elseif fr<f(n)
  sim(end,:)=xr;f(end)=fr;
 elseif fr<f(end)
  xc=min(max(1.5*center-.5*sim(end,:),lower),upper);fc=fun(xc);
  if fc<=fr,sim(end,:)=xc;f(end)=fc;else,shrink=true;end
 else
  xc=min(max(.5*center+.5*sim(end,:),lower),upper);fc=fun(xc);
  if fc<f(end),sim(end,:)=xc;f(end)=fc;else,shrink=true;end
 end
 if shrink
  for j=2:n+1,sim(j,:)=min(max(sim(1,:)+.5*(sim(j,:)-sim(1,:)),lower),upper);f(j)=fun(sim(j,:));end
 end
 [f,ix]=sort(f);sim=sim(ix,:);iterations=iterations+1;
end
q=sim(1,:);best=f(1);
end
