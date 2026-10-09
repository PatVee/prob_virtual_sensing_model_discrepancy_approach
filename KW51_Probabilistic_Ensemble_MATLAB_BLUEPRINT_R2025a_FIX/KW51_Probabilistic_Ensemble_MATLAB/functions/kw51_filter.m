function y=kw51_filter(x,cutoff,fs)
% Fourth-order Butterworth, zero phase, 15-sample odd-reflection padding.
% Matches the preprocessing used in the Python study; no fitted values read.
[b,a]=butter(4,cutoff/(fs/2));b=b/a(1);a=a/a(1);n=4;edge=15;
assert(size(x,1)>edge,'Record too short for filtering.');
A=[-a(2:end).',[eye(n-1);zeros(1,n-1)]];
zi=(eye(n)-A)\(b(2:end).'-a(2:end).'*b(1));
y=zeros(size(x));
for j=1:size(x,2)
 z=[2*x(1,j)-x(edge+1:-1:2,j);x(:,j);2*x(end,j)-x(end-1:-1:end-edge,j)];
 z=filter(b,a,z,zi*z(1));z=flipud(z);z=filter(b,a,z,zi*z(1));z=flipud(z);
 y(:,j)=z(edge+1:end-edge);
end
end
