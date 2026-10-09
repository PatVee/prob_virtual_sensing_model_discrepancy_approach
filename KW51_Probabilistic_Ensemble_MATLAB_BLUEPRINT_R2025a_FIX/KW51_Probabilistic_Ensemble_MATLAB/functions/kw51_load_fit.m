function [loads,score]=kw51_load_fit(H,y,W)
n=size(H,3);T=size(H,1);R=size(H,2);X=zeros(T*R,n);
for j=1:n,z=H(:,:,j)*W.';X(:,j)=z(:);end
z=y*W.';z=z(:);J=eye(n)-ones(n)/n;
reg=J/.30^2+ones(n)/(n*n*.45^2);
A=X.'*X+reg;b=X.'*z+ones(n,1)*1.4/(n*.45^2);
[loads,~]=kw51_boxqp(A,b,.6,2.2);
err=X*loads-z;
score=(err.'*err+loads.'*reg*loads-2*sum(loads)*1.4/(n*.45^2)+1.4^2/.45^2)/numel(z);
end
