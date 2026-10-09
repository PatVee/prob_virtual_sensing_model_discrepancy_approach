function [x,kkt]=kw51_boxqp(A,b,lo,hi)
% Strictly convex, box-constrained quadratic minimisation without a toolbox.
% Primal active-set method: minimise 0.5*x'*A*x-b'*x.
n=numel(b);x=min(max(A\b,lo),hi);tol=1e-9;
bound=zeros(n,1);bound(x<=lo+tol)=-1;bound(x>=hi-tol)=1;
for it=1:20*n+100
 free=find(bound==0);fixed=find(bound~=0);z=x;
 if ~isempty(free),z(free)=A(free,free)\(b(free)-A(free,fixed)*x(fixed));end
 p=z-x;
 if max(abs(p))<1e-10
  g=A*x-b;violation=zeros(n,1);
  violation(bound==-1)=max(0,-g(bound==-1));violation(bound==1)=max(0,g(bound==1));
  [bad,j]=max(violation);
  if bad<1e-7,break;end
  bound(j)=0;
 else
  step=1;hit=0;hitbound=0;
  for j=free(:).'
   if p(j)>0,v=(hi-x(j))/p(j);s=1;elseif p(j)<0,v=(lo-x(j))/p(j);s=-1;else,continue;end
   if v<step,step=max(0,v);hit=j;hitbound=s;end
  end
  x=min(max(x+step*p,lo),hi);
  if hit>0,bound(hit)=hitbound;end
 end
end
g=A*x-b;pg=g;pg(x<=lo+1e-8)=min(g(x<=lo+1e-8),0);pg(x>=hi-1e-8)=max(g(x>=hi-1e-8),0);kkt=max(abs(pg));
assert(kkt<1e-4,'KW51:LoadSolve','Bounded load solve failed its optimality check.');
end
