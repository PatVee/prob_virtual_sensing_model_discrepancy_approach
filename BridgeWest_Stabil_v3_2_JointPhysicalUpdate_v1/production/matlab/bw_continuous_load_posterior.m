function out=bw_continuous_load_posterior(Bref_microstrain_per_N,yref,discrepancy,c,varargin)
%BW_CONTINUOUS_LOAD_POSTERIOR Bounded 8-bogie posterior from reference strain only.
% Load coordinates are kN. The likelihood covariance is supplied by the
% shaker-informed discrepancy object and is common to all load states.
ip=inputParser;
addParameter(ip,'Seed',c.gibbs_seed,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'Burn',c.gibbs_burn_final,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'Keep',c.gibbs_keep_final,@(x)isnumeric(x)&&isscalar(x));
parse(ip,varargin{:});opt=ip.Results;
T=size(Bref_microstrain_per_N,1);R=size(Bref_microstrain_per_N,2);nb=size(Bref_microstrain_per_N,3);
assert(nb==8 && size(yref,1)==T && size(yref,2)==R);
H=zeros(nb);h=zeros(nb,1);ti=1:c.likelihood_time_stride:T;
for jt=1:numel(ti)
    it=ti(jt);
    B=reshape(Bref_microstrain_per_N(it,:,:),R,nb)*1000; % microstrain/kN
    S=discrepancy.reference_precision(:,:,it);
    d=yref(it,:).'-discrepancy.mean{1}(it,discrepancy.reference_indices).';
    H=H+c.likelihood_power*(B.'*S*B);
    h=h+c.likelihood_power*(B.'*S*d);
end
lo=repmat(double(c.load_lower_kN),nb,1);hi=repmat(double(c.load_upper_kN),nb,1);
g=bw_bounded_gaussian_moments(H,h,lo,hi,'Seed',opt.Seed,'Burn',opt.Burn,'Keep',opt.Keep);
out=g;out.time_indices=ti;out.condition_number=cond(g.precision);out.lower_kN=lo;out.upper_kN=hi;
end
