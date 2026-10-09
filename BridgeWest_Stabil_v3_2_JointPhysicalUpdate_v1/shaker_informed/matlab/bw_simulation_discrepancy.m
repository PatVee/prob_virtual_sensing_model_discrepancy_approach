function discrepancy=bw_simulation_discrepancy(delta,refs,c)
% Compute original simulation discrepancy moments without measured responses.
M=size(delta{1},1); T=size(delta{1},2); R=numel(refs);
mu=cell(1,3); maps=cell(1,3); cv=cell(1,3); invR=zeros(R,R,T);
for k=1:3
 O=size(delta{k},3); mu{k}=reshape(mean(delta{k},1),T,O);
 maps{k}=zeros(O,R,T); cv{k}=zeros(T,O);
end
for it=1:T
 D=reshape(delta{1}(:,it,refs),M,R); D=D-mean(D,1);
 C=D.'*D/(M-1); C=(C+C.')/2;
 S=C+eye(R)*(c.cov_ridge*max(trace(C)/R,1e-12)+1e-12+c.reference_noise_sigma_microstrain^2);
 invR(:,:,it)=pinv(S,1e-15*norm(S,2));
 for k=1:3
  O=size(delta{k},3); A=reshape(delta{k}(:,it,:),M,O); A=A-mean(A,1);
  cross=A.'*D/(M-1); maps{k}(:,:,it)=cross*invR(:,:,it);
  variance=sum(A.^2,1)/(M-1)-sum(maps{k}(:,:,it).*cross,2).';
  if k==1, noise=c.target_noise_sigma_microstrain; else, noise=c.displacement_noise_sigma_mm; end
  cv{k}(it,:)=max(variance,0)+noise^2;
 end
end
discrepancy=struct('mean',{mu},'conditional_gain',{maps},'conditional_variance',{cv},'reference_precision',invR,'reference_indices',refs,'paired_draw_count',M);
end
