function z=bw_truncnorm_scalar(a,b)
%BW_TRUNCNORM_SCALAR Draw N(0,1) truncated to [a,b], without Statistics Toolbox.
assert(isfinite(a)||isinf(a)); assert(isfinite(b)||isinf(b)); assert(a<b);
% Symmetry keeps the rejection sampler in the positive tail.
if b < -2
    z=-bw_truncnorm_scalar(-b,-a); return;
end
% Narrow finite intervals: uniform proposal with exact density rejection.
if isfinite(a) && isfinite(b) && (b-a)<1.0
    if a<=0 && b>=0, z0=0; else, z0=min(abs([a b])); end
    logM=-0.5*z0^2;
    for k=1:100000
        q=a+(b-a)*rand;
        if log(rand) <= -0.5*q^2-logM, z=q; return; end
    end
    error('BridgeWest:TruncNormal','Narrow-interval rejection sampler did not accept.');
end
% Robert exponential proposal for a positive normal tail.
if a > 2
    alpha=(a+sqrt(a*a+4))/2;
    for k=1:100000
        q=a-log(rand)/alpha;
        if q>b, continue; end
        if log(rand) <= -0.5*(q-alpha)^2, z=q; return; end
    end
    error('BridgeWest:TruncNormal','Positive-tail rejection sampler did not accept.');
end
% Central interval inverse CDF. erfc/erfcinv are base MATLAB functions.
if isinf(a), Fa=0; else, Fa=0.5*erfc(-a/sqrt(2)); end
if isinf(b), Fb=1; else, Fb=0.5*erfc(-b/sqrt(2)); end
p=Fa+(Fb-Fa)*rand;
p=min(max(p,realmin('double')),1-eps('double'));
z=-sqrt(2)*erfcinv(2*p);
% Roundoff guard.
z=min(max(z,a),b);
end
