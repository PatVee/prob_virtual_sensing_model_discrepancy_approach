function Q=kwpe_condition(nom,delta,y,refs,c,t)
% Same Gaussian cross-covariance and mixture moments as bw_v2_predict.
N=size(nom,1);T=size(nom,2);O=size(nom,3);M=size(delta,1);R=numel(refs);
assert(isequal(size(y),[T R])&&all(isfinite(y(:))),'Invalid reference observations.');
mu=reshape(mean(delta,1),T,O);dc=delta-reshape(mu,1,T,O);
gain=zeros(O,R,T);conditionalVariance=zeros(T,O);inverse=zeros(R,R,T);refNoise=zeros(T,1);
outputNoise=[repmat(c.target_noise_sigma_microstrain,1,5),repmat(c.displacement_noise_sigma_mm,1,2)];
assert(O==7,'KW51 output order requires five strains and two displacements.');
for it=1:T
    D=reshape(dc(:,it,:),M,O);C=(D.'*D)/(M-1);C=(C+C.')/2;
    ridge=c.cov_ridge*max(trace(C(refs,refs))/R,1e-12)+1e-12;
    refNoise(it)=c.reference_noise_sigma_microstrain^2+ridge;
    S=C(refs,refs)+eye(R)*refNoise(it);inverse(:,:,it)=S\eye(R);
    gain(:,:,it)=C(:,refs)*inverse(:,:,it);
    V=C-gain(:,:,it)*C(refs,:);
    conditionalVariance(it,:)=max(diag(V).',0)+outputNoise.^2;
end
ll=zeros(N,1);likelihoodIndices=1:c.likelihood_time_stride:T;
for it=likelihoodIndices
    innovation=y(it,:)-reshape(nom(:,it,refs),N,R)-mu(it,refs);
    ll=ll-.5*c.likelihood_power*sum((innovation*inverse(:,:,it)).*innovation,2);
end
weights=exp(ll-max(ll));weights=weights/sum(weights);
mn=zeros(T,O);variance=zeros(T,O);
for it=1:T
    prediction=reshape(nom(:,it,:),N,O);
    innovation=y(it,:)-prediction(:,refs)-mu(it,refs);
    conditioned=prediction+mu(it,:)+innovation*gain(:,:,it).';
    mn(it,:)=weights.'*conditioned;
    variance(it,:)=max(weights.'*(conditioned.^2)-mn(it,:).^2+conditionalVariance(it,:),0);
end
S=c.posterior_draw_count;edges=cumsum(weights);edges(end)=1;
u=rand(S,1);selected=zeros(S,1);
for s=1:S,selected(s)=find(edges>=u(s),1,'first');end
Z=randn(S,M)/sqrt(M-1);
residual=reshape(Z*reshape(dc,M,T*O),S,T,O);
trajectories=zeros(S,T,O);
for it=1:T
    prediction=reshape(nom(selected,it,:),S,O);
    innovation=y(it,:)-prediction(:,refs)-mu(it,refs);
    a=prediction+mu(it,:)+innovation*gain(:,:,it).';
    residualAtTime=reshape(residual(:,it,:),S,O);
    referenceNoise=randn(S,R)*sqrt(refNoise(it));
    r=residualAtTime-(residualAtTime(:,refs)+referenceNoise)*gain(:,:,it).';
    r=r+randn(S,O).*outputNoise;
    trajectories(:,it,:)=reshape(a+r,S,1,O);
end
Q=struct('mean',mn,'structural_sd',sqrt(variance),'weights',weights, ...
    'posterior_ess',1/sum(weights.^2),'gain',gain,'conditional_variance',conditionalVariance, ...
    'prior_mean',reshape(mean(nom,1),T,O)+mu,'trajectories',trajectories, ...
    'likelihood_indices',likelihoodIndices);
end
