function nom=kwmo_nominal(B,E,F,loads)
% Fixed unscaled quasi-static FE response; candidate variability is load only.
H=kw51_basis(B,E,E.time_s(:),F.kinematics(:).');
T=size(H,1);O=size(H,2);na=size(H,3);
assert(size(loads,2)==na&&all(isfinite(loads(:))),'Invalid axle-load ensemble.');
nom=zeros(size(loads,1),T,O);
for o=1:O
 h=reshape(H(:,o,:),T,na);nom(:,:,o)=loads*h.';
end
end
