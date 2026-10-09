function [p,discrepancy]=bw_predict_shaker(nom,delta,deltaWeights,y,refs,target,c)
%BW_PREDICT_SHAKER Reference-conditioned prediction with weighted shaker discrepancy.
% nom: cell with N x T x outputs nominal candidates.
% delta: cell with M x T x outputs paired Stabil shaker discrepancies.
% deltaWeights: M x 1 discrete shaker-state weights.
% y: T x R measured reference strains only.
N=size(nom{1},1);T=size(nom{1},2);R=numel(refs);
assert(size(y,2)==R && all(isfinite(y(:))),'Reference observations must be finite.');
discrepancy=bw_simulation_discrepancy_weighted(delta,deltaWeights,refs,c);
mu=discrepancy.mean;maps=discrepancy.conditional_gain;cv=discrepancy.conditional_variance;invR=discrepancy.reference_precision;
ll=zeros(N,1);ti=1:c.likelihood_time_stride:T;
for it=ti
    inn=y(it,:)-double(reshape(nom{1}(:,it,refs),N,R))-mu{1}(it,refs);
    ll=ll-.5*c.likelihood_power*sum((inn*invR(:,:,it)).*inn,2);
end
w=exp(ll-max(ll));w=w/sum(w);p.posterior_ess=1/sum(w.^2);
p.posterior_weights=w;p.likelihood_time_indices=ti;
keys={'strain','disp_BA','disp_BB'};
for k=1:3
    O=size(nom{k},3);mn=zeros(T,O);sd=mn;
    for it=1:T
        inn=y(it,:)-double(reshape(nom{1}(:,it,refs),N,R))-mu{1}(it,refs);
        a=double(reshape(nom{k}(:,it,:),N,O))+mu{k}(it,:)+inn*maps{k}(:,:,it).';
        mn(it,:)=w.'*a;sd(it,:)=sqrt(max(w.'*(a.^2)-mn(it,:).^2+cv{k}(it,:),0));
    end
    if k==1,mn=mn(:,target);sd=sd(:,target);end
    p.([keys{k} '_mean'])=mn;p.([keys{k} '_sd'])=sd;
    p.([keys{k} '_lo'])=mn-c.interval_z*sd;p.([keys{k} '_hi'])=mn+c.interval_z*sd;
end
end
