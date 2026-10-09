function discrepancy=bw_simulation_discrepancy_weighted(delta,weights,refs,c)
%BW_SIMULATION_DISCREPANCY_WEIGHTED Weighted paired-simulation discrepancy.
% weights define the discrete shaker-state mixture and must sum to one.
M=size(delta{1},1);T=size(delta{1},2);R=numel(refs);
w=double(weights(:));assert(numel(w)==M && all(w>=0),'Invalid discrepancy weights.');
w=w/sum(w);
mu=cell(1,3);maps=cell(1,3);cv=cell(1,3);invR=zeros(R,R,T);
for k=1:3
    O=size(delta{k},3);A=reshape(delta{k},M,[]);m=w.'*A;mu{k}=reshape(m,T,O);
    maps{k}=zeros(O,R,T);cv{k}=zeros(T,O);
end
for it=1:T
    D=reshape(delta{1}(:,it,refs),M,R);D=D-mu{1}(it,refs);
    C=D.'*(D.*w);C=(C+C.')/2;
    ridge=c.cov_ridge*max(trace(C)/R,1e-12)+1e-12+c.reference_noise_sigma_microstrain^2;
    S=C+eye(R)*ridge;invR(:,:,it)=pinv(S,1e-15*norm(S,2));
    for k=1:3
        O=size(delta{k},3);A=reshape(delta{k}(:,it,:),M,O);A=A-mu{k}(it,:);
        cross=A.'*(D.*w);maps{k}(:,:,it)=cross*invR(:,:,it);
        variance=sum((A.^2).*w,1)-sum(maps{k}(:,:,it).*cross,2).';
        if k==1,noise=c.target_noise_sigma_microstrain;else,noise=c.displacement_noise_sigma_mm;end
        cv{k}(it,:)=max(variance,0)+noise^2;
    end
end
discrepancy=struct('mean',{mu},'conditional_gain',{maps},'conditional_variance',{cv}, ...
    'reference_precision',invR,'reference_indices',refs,'paired_draw_count',M, ...
    'state_weights',w,'effective_state_count',1/sum(w.^2));
end
