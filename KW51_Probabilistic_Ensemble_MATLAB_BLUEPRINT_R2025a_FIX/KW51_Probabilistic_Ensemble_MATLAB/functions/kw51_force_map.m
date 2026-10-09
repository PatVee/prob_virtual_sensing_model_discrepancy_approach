function [G,W,K]=kw51_force_map(B,refs,forceScale)
G=full(B.G)*1e5;gr=G(refs,:);
C=.5^2*eye(numel(refs))+forceScale^2*(gr*gr.');
[V,D]=eig((C+C.')/2);W=V*diag(1./sqrt(diag(D)))*V.';
K=forceScale^2*gr.'/C;
end
