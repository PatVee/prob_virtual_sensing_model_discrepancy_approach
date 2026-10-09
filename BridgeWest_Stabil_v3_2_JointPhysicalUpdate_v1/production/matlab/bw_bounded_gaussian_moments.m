function out=bw_bounded_gaussian_moments(H,h,lo,hi,varargin)
%BW_BOUNDED_GAUSSIAN_MOMENTS Gibbs moments for exp(-.5*x'Hx+h'x) on a box.
% No Statistics/Optimization Toolbox is required.
ip=inputParser;
addParameter(ip,'Seed',20260924,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'Burn',1000,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(ip,'Keep',8000,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(ip,'Thin',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
parse(ip,varargin{:}); opt=ip.Results;
H=(double(H)+double(H).')/2; h=double(h(:)); lo=double(lo(:)); hi=double(hi(:));
n=numel(h); assert(all(size(H)==[n n])&&numel(lo)==n&&numel(hi)==n&&all(lo<hi));
ev=eig(H); target=1e-12*max(max(ev),1); ridge=max(0,target-min(ev));
if ridge>0, H=H+ridge*eye(n); ev=eig(H); end
m=H\h;
% Deterministic box-QP coordinate descent gives a stable Gibbs initial state.
x=min(max(m,lo),hi);
for sweep=1:10000
    old=x;
    for i=1:n
        hii=H(i,i);
        cm=m(i)-(H(i,:)*(x-m)-hii*(x(i)-m(i)))/hii;
        x(i)=min(max(cm,lo(i)),hi(i));
    end
    if max(abs(x-old))<1e-11, break; end
end
map=x;
oldrng=rng; cleaner=onCleanup(@()rng(oldrng)); rng(opt.Seed,'twister');
S=zeros(opt.Keep,n); nkeep=0; total=opt.Burn+opt.Keep*opt.Thin;
for sweep=1:total
    for i=1:n
        hii=H(i,i); sd=sqrt(1/hii);
        cm=m(i)-(H(i,:)*(x-m)-hii*(x(i)-m(i)))/hii;
        a=(lo(i)-cm)/sd; b=(hi(i)-cm)/sd;
        x(i)=cm+sd*bw_truncnorm_scalar(a,b);
    end
    if sweep>opt.Burn && mod(sweep-opt.Burn-1,opt.Thin)==0
        nkeep=nkeep+1; S(nkeep,:)=x.';
    end
end
clear cleaner
assert(nkeep==opt.Keep);
mu=mean(S,1).'; C=cov(S,0); if n==1,C=reshape(C,1,1);end
out=struct('mean',mu,'covariance',C,'samples',S,'unconstrained_mean',m,'box_map',map, ...
    'precision',H,'precision_eigenvalues',ev,'ridge',ridge,'seed',opt.Seed,'burn',opt.Burn,'keep',opt.Keep,'thin',opt.Thin);
end
