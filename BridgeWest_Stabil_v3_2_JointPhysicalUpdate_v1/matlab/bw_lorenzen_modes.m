function fe=bw_lorenzen_modes(fe)
% No mass is added to massless shell drilling rotations.
% Shift-invert operator supports semidefinite physical mass matrices.
f=fe.free;K=fe.K(f,f);M=fe.M(f,f);nm=fe.c.number_of_modes;
sigma=-1e-4;solver=decomposition(K-sigma*M,'lu');
opts=struct('isreal',true,'issym',false,'tol',2e-9,'maxit',2000);
opts.v0=sin((1:numel(f))');
[V,D,flag]=eigs(@(x)solver\(M*x),numel(f),nm,'largestabs',opts);
assert(flag==0,'Modal eigensolver did not converge.');
values=sigma+1./diag(D);
assert(max(abs(imag(values)))<1e-6*max(abs(values)),'Non-real eigenvalues.');
[values,ix]=sort(real(values));V=real(V(:,ix));
assert(all(values>0),'Nonpositive structural modes: check restraints.');
V=V./sqrt(sum(V.*(M*V),1));
P=zeros(size(fe.K,1),nm);P(f,:)=V;
fe.lam=values;fe.freq=sqrt(values)/(2*pi);fe.P=P;
fe.residual=sqrt(sum((K*V-(M*V).*values').^2,1))./(sqrt(sum((K*V).^2,1))+sqrt(sum(((M*V).*values').^2,1)));
fe.mass_orthogonality_error=norm(V'*M*V-eye(nm),'fro');
end
