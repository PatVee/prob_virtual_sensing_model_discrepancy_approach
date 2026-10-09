function D=kwmo_decompose(nom,delta,y,Q,c)
% Exact pointwise mixture decomposition using the SAME posterior weights.
N=size(nom,1);T=size(nom,2);O=size(nom,3);w=Q.weights(:);
mu=reshape(mean(delta,1),T,O);
noise=[repmat(c.target_noise_sigma_microstrain^2,1,5),repmat(c.displacement_noise_sigma_mm^2,1,2)];
fields={'fe_mean','correction_mean','fe_variance','correction_between_variance','twice_fe_correction_covariance'};
for j=1:numel(fields),D.(fields{j})=zeros(T,O);end
for it=1:T
 a=reshape(nom(:,it,:),N,O);
 b=mu(it,:)+(y(it,:)-a(:,1:3)-mu(it,1:3))*Q.gain(:,:,it).';
 am=w.'*a;bm=w.'*b;ac=a-am;bc=b-bm;
 D.fe_mean(it,:)=am;D.correction_mean(it,:)=bm;
 D.fe_variance(it,:)=w.'*(ac.^2);
 D.correction_between_variance(it,:)=w.'*(bc.^2);
 D.twice_fe_correction_covariance(it,:)=2*w.'*(ac.*bc);
end
D.conditional_discrepancy_variance=max(Q.conditional_variance-noise,0);
D.output_noise_variance=repmat(noise,T,1);
D.total_variance=D.fe_variance+D.correction_between_variance+ ...
 D.twice_fe_correction_covariance+D.conditional_discrepancy_variance+D.output_noise_variance;
D.corrected_mean=D.fe_mean+D.correction_mean;
D.mean_identity_error=max(abs(D.corrected_mean(:)-Q.mean(:)));
D.variance_identity_error=max(abs(D.total_variance(:)-Q.structural_sd(:).^2));
assert(D.mean_identity_error<1e-8*max(1,max(abs(Q.mean(:)))));
assert(D.variance_identity_error<1e-8*max(1,max(Q.structural_sd(:).^2)));
D.note='FE variance comes from candidate axle loads at one fixed FE baseline. Covariance can be negative; components are not independent percentages. Pointwise unfiltered decomposition.';
end
