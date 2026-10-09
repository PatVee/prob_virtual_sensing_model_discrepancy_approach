function delta=bw_project_discrepancy_basis(deltaBasis,loadsN)
%BW_PROJECT_DISCREPANCY_BASIS Project unit-bogie discrepancy bases to loads.
% deltaBasis{k}: M x T x O x 8. loadsN may be 1 x 8 (same load for every
% shaker state) or M x 8 (state-specific quadrature/loading design).
M=size(deltaBasis{1},1);loadsN=double(loadsN);
assert(size(loadsN,2)==8,'loadsN must have 8 bogie loads.');
if size(loadsN,1)==1,loadsN=repmat(loadsN,M,1);end
assert(size(loadsN,1)==M,'loadsN must be 1 x 8 or M x 8.');
delta=cell(1,3);
for k=1:3
    T=size(deltaBasis{k},2);O=size(deltaBasis{k},3);delta{k}=zeros(M,T,O);
    for m=1:M
        D=squeeze(deltaBasis{k}(m,:,:,:)); % T x O x 8
        delta{k}(m,:,:)=reshape(sum(D.*reshape(loadsN(m,:),1,1,8),3),1,T,O);
    end
end
end
