function F=bw_stabil_force_centerline(fe,t,v,c,direction)
%BW_STABIL_FORCE_CENTERLINE Eight unit-bogie moving-load bases on the deck.
% Each bogie contains two axles. Each axle carries one half of a unit bogie
% resultant and is applied at the track centreline. No unverified rail spacing
% or wheel-load imbalance is introduced in this discrepancy stage.
if nargin<5 || isempty(direction),direction=1;end
assert(isscalar(direction)&&any(direction==[-1 1]),'direction must be +1 or -1.');
ax=double(c.axle_offsets_m(:));assert(numel(ax)==16,'Expected 16 axle offsets.');
F=zeros(size(fe.K,1),8);a=fe.c.bearing_x_m(1);b=fe.c.bearing_x_m(2);
if direction>0,xlead=a+v*t;else,xlead=b-v*t;end
ytrack=mean(fe.c.girder_y_m);
for ia=1:numel(ax)
    x=xlead-direction*ax(ia);
    if x<a-1e-12 || x>b+1e-12,continue;end
    bogie=ceil(ia/2);
    F(:,bogie)=addDeckPoint(F(:,bogie),x,ytrack,-0.5);
end
    function q=addDeckPoint(q,xp,yp,forceN)
        candidates=find(fe.shells(:,6)==4);
        pick=0;best=inf;
        for ee=candidates(:)'
            ids=fe.shells(ee,1:4);p=fe.nodes(ids,:);
            xmin=min(p(:,1));xmax=max(p(:,1));ymin=min(p(:,2));ymax=max(p(:,2));
            if xp>=xmin-1e-10 && xp<=xmax+1e-10 && yp>=ymin-1e-10 && yp<=ymax+1e-10
                d=(xp-mean(p(:,1)))^2+(yp-mean(p(:,2)))^2;
                if d<best,best=d;pick=ee;end
            end
        end
        assert(pick>0,'Moving load lies outside deck shell mesh.');
        ids=fe.shells(pick,1:4);p=fe.nodes(ids,:);
        xmin=min(p(:,1));xmax=max(p(:,1));ymin=min(p(:,2));ymax=max(p(:,2));
        xi=2*(xp-xmin)/(xmax-xmin)-1;eta=2*(yp-ymin)/(ymax-ymin)-1;
        N=.25*[(1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];
        dof=6*(ids-1)+3;
        q(dof(:))=q(dof(:))+forceN*N(:);
    end
end
