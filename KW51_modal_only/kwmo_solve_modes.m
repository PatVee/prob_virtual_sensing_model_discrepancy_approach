function [phi,omega,checks]=kwmo_solve_modes(K,M,nModes)
% Diagonal congruence scaling and inverse-subspace refinement for KW51 FE.
K=(K+K')/2;M=(M+M')/2;n=size(K,1);
assert(all(diag(K)>0),'Nonpositive stiffness diagonal.');
s=1./sqrt(diag(K));D=spdiags(s,0,n,n);Ks=D*K*D;Ms=D*M*D;
opts.tol=1e-11;opts.maxit=5000;opts.p=min(n,max(4*nModes,100));
opts.v0=ones(n,1);
[V,~,flag]=eigs(Ks,Ms,nModes,0,opts);
assert(flag==0,'Modal eigensolver did not converge.');phi=D*V;
for j=1:2
 W=D*(Ks\(D*(M*phi)));
 Gram=full(W'*M*W);[R,p]=chol((Gram+Gram')/2);
 assert(p==0,'Modal subspace has a nonpositive mass Gram matrix.');
 V=W/R;Kr=full(V'*K*V);[U,L]=eig((Kr+Kr')/2);
 [lambda,ix]=sort(real(diag(L)));phi=V*U(:,ix);
end
assert(isreal(phi)&&all(lambda>0),'Nonpositive/complex baseline modes.');
omega=sqrt(lambda);res=K*phi-M*phi*diag(lambda);
checks.mass_orthogonality_error=norm(phi'*M*phi-eye(nModes),'fro');
checks.eigen_residual=norm(res,'fro')/max(norm(K*phi,'fro'),eps);
checks.eigen_backward_error=max(sqrt(sum(abs(res).^2,1))./ ...
 ((norm(K,'fro')+lambda.'*norm(M,'fro')).*sqrt(sum(abs(phi).^2,1))));
assert(checks.mass_orthogonality_error<1e-5,'Mass normalization failed.');
assert(checks.eigen_residual<1e-5,'Eigenproblem residual too large: %.3g.',checks.eigen_residual);
assert(checks.eigen_backward_error<1e-12,'Eigenproblem backward error too large.');
end
