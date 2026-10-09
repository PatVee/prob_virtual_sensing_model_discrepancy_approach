function [nom,B0,meta]=bw_stabil_nominal_candidates(fe,v,tobs,loadsN,c,direction)
%BW_STABIL_NOMINAL_CANDIDATES Build the baseline Stabil moving-load library.
% loadsN is N x 8 bogie loads in N. Structural variability is represented by
% the shaker discrepancy ensemble, so no additional arbitrary global stiffness
% grid is applied here by default.
if nargin<6 || isempty(direction),direction=1;end
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<5 || isempty(c),c=jsondecode(fileread(fullfile(root,'discrepancy_config.json')));end
ens=bw_shaker_ensemble(root);[B0,meta]=bw_stabil_basis_shaker(fe,v,ens.baseline,tobs,c,direction);
loadsN=double(loadsN);assert(size(loadsN,2)==8 && all(isfinite(loadsN(:))),'loadsN must be finite N x 8.');
N=size(loadsN,1);T=numel(tobs);nom={zeros(N,T,7,'single'),zeros(N,T,41,'single'),zeros(N,T,41,'single')};
for n=1:N
    l3=reshape(loadsN(n,:),1,1,8);
    for k=1:3
        y=sum(B0{k}.*l3,3);nom{k}(n,:,:)=single(reshape(y,1,size(y,1),size(y,2)));
    end
end
end
