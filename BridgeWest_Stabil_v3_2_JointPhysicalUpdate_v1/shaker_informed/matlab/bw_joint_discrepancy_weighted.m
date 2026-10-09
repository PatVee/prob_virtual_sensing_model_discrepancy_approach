function D=bw_joint_discrepancy_weighted(components,loadsN,refs,c)
%BW_JOINT_DISCREPANCY_WEIGHTED Independent multi-source discrepancy moments.
%
% components is a struct array with fields:
%   name        descriptive label
%   delta_basis 1x3 cell; each S x T x O x 8
%   weights     Sx1 non-negative weights
%
% Each component is projected through the same 1x8 load vector, then the
% component means and covariances are added under an independence assumption.
% This keeps shaker dynamics, load transfer, common compliance and strain
% recovery as distinct uncertainty sources without constructing their full
% Cartesian product.
assert(~isempty(components),'At least one discrepancy component is required.');
loadsN=double(loadsN(:).');assert(numel(loadsN)==8,'loadsN must contain eight bogie loads.');
R=numel(refs);T=size(components(1).delta_basis{1},2);
mu={zeros(T,size(components(1).delta_basis{1},3)), ...
    zeros(T,size(components(1).delta_basis{2},3)), ...
    zeros(T,size(components(1).delta_basis{3},3))};
Cen=cell(numel(components),3);W=cell(numel(components),1);
for ic=1:numel(components)
    delta=bw_project_discrepancy_basis(components(ic).delta_basis,loadsN);
    w=double(components(ic).weights(:));assert(numel(w)==size(delta{1},1)&&all(w>=0),'Invalid component weights.');
    w=w/sum(w);W{ic}=w;
    for k=1:3
        S=size(delta{k},1);O=size(delta{k},3);
        A=reshape(delta{k},S,[]);mk=w.'*A;mk=reshape(mk,T,O);
        mu{k}=mu{k}+mk;
        Cen{ic,k}=delta{k}-reshape(mk,1,T,O);
    end
end
maps=cell(1,3);cv=cell(1,3);invR=zeros(R,R,T);CrrAll=zeros(R,R,T);
for k=1:3
    maps{k}=zeros(size(mu{k},2),R,T);cv{k}=zeros(T,size(mu{k},2));
end
for it=1:T
    Crr=zeros(R);
    for ic=1:numel(components)
        Dr=reshape(Cen{ic,1}(:,it,refs),[],R);w=W{ic};
        Crr=Crr+Dr.'*(Dr.*w);
    end
    Crr=(Crr+Crr.')/2;CrrAll(:,:,it)=Crr;
    ridge=c.cov_ridge*max(trace(Crr)/R,1e-12)+1e-12+c.reference_noise_sigma_microstrain^2;
    S=Crr+eye(R)*ridge;invR(:,:,it)=pinv(S,1e-15*norm(S,2));
    for k=1:3
        O=size(mu{k},2);cross=zeros(O,R);variance=zeros(1,O);
        for ic=1:numel(components)
            A=reshape(Cen{ic,k}(:,it,:),[],O);Dr=reshape(Cen{ic,1}(:,it,refs),[],R);w=W{ic};
            cross=cross+A.'*(Dr.*w);variance=variance+sum((A.^2).*w,1);
        end
        maps{k}(:,:,it)=cross*invR(:,:,it);
        variance=variance-sum(maps{k}(:,:,it).*cross,2).';
        if k==1,noise=c.target_noise_sigma_microstrain;else,noise=c.displacement_noise_sigma_mm;end
        cv{k}(it,:)=max(variance,0)+noise^2;
    end
end
D=struct('mean',{mu},'conditional_gain',{maps},'conditional_variance',{cv}, ...
    'reference_precision',invR,'reference_covariance',CrrAll,'reference_indices',refs, ...
    'component_names',{cellstr(string({components.name}))});
end
