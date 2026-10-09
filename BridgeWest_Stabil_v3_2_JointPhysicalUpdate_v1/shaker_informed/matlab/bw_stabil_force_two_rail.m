function F=bw_stabil_force_two_rail(fe,t,v,c,direction)
%BW_STABIL_FORCE_TWO_RAIL Explicit two-wheel/rail moving-load transfer.
%
% For each bogie basis P_b = 1 N:
%   two axles              -> P_b/2 per axle
%   two wheel/rail lines   -> P_b/4 per wheel
%   rail/ballast footprint -> distributed vertical nodal forces on deck
%
% The footprint is integrated over tributary rectangles of the deck-shell
% nodes. When conserve_wheel_resultant=true (default), the nodal weights for
% each wheel are renormalised so their sum is one while the wheel centre is
% within the original bearing-to-bearing moving-load interval. Consequently
% the new operator changes only the spatial force transfer; it does not alter
% the total applied bogie resultant or passage timing.
%
% Required/optional fields in c:
%   axle_offsets_m
%   rail_gauge_m                              default 1.435
%   rail_spread_longitudinal_half_length_m    default 2.4
%   rail_spread_transverse_half_width_m       default 0.6
%   conserve_wheel_resultant                  default true
%
% NOTE: the rail/ballast footprint dimensions are modelling assumptions and
% are deliberately exposed in discrepancy_config.json for sensitivity work.
if nargin<5 || isempty(direction),direction=1;end
assert(isscalar(direction)&&any(direction==[-1 1]),'direction must be +1 or -1.');
ax=double(c.axle_offsets_m(:));assert(numel(ax)==16,'Expected 16 axle offsets.');

gauge=getcfg(c,'rail_gauge_m',1.435);
hx=getcfg(c,'rail_spread_longitudinal_half_length_m',2.4);
hy=getcfg(c,'rail_spread_transverse_half_width_m',0.6);
conserve=getcfg(c,'conserve_wheel_resultant',true);
assert(gauge>0 && hx>0 && hy>0,'Two-rail spread dimensions must be positive.');

F=zeros(size(fe.K,1),8);
a=fe.c.bearing_x_m(1);b=fe.c.bearing_x_m(2);
if direction>0,xlead=a+v*t;else,xlead=b-v*t;end
ytrack=mean(fe.c.girder_y_m);
yrail=[ytrack-gauge/2,ytrack+gauge/2];

G=deckGrid(fe);
for ia=1:numel(ax)
    x=xlead-direction*ax(ia);
    bogie=ceil(ia/2);

    % Preserve the original moving-load timing when conserving the wheel
    % resultant: the axle becomes active at the same bearing coordinate as
    % the former centreline point load.
    if conserve && (x<a-1e-12 || x>b+1e-12),continue;end

    for ir=1:2
        w=footprintWeights(G,x,yrail(ir),hx,hy,a,b,conserve);
        if isempty(w),continue;end
        dof=6*(G.ids-1)+3;
        % Unit bogie -> 1/2 per axle -> 1/4 per wheel.
        F(dof,bogie)=F(dof,bogie)-0.25*w;
    end
end
end

function v=getcfg(c,name,default)
if isfield(c,name) && ~isempty(c.(name)),v=c.(name);else,v=default;end
end

function G=deckGrid(fe)
% Build nodal tributary rectangles for deck-shell group 4. A persistent cache
% avoids rebuilding the grid at every Newmark time step.
persistent key cache
newkey=[size(fe.nodes,1),size(fe.shells,1), ...
    min(fe.nodes(:,1)),max(fe.nodes(:,1)),min(fe.nodes(:,2)),max(fe.nodes(:,2))];
if ~isempty(key) && isequal(key,newkey)
    G=cache;return
end

deck=find(fe.shells(:,6)==4);
ids=unique(fe.shells(deck,1:4));ids=ids(:);
x=fe.nodes(ids,1);y=fe.nodes(ids,2);
ux=unique(round(x,10));uy=unique(round(y,10));
xb=zeros(numel(ux),2);yb=zeros(numel(uy),2);
for i=1:numel(ux)
    if i==1,xb(i,1)=ux(i);else,xb(i,1)=(ux(i-1)+ux(i))/2;end
    if i==numel(ux),xb(i,2)=ux(i);else,xb(i,2)=(ux(i)+ux(i+1))/2;end
end
for i=1:numel(uy)
    if i==1,yb(i,1)=uy(i);else,yb(i,1)=(uy(i-1)+uy(i))/2;end
    if i==numel(uy),yb(i,2)=uy(i);else,yb(i,2)=(uy(i)+uy(i+1))/2;end
end
xL=zeros(size(x));xR=xL;yB=xL;yT=xL;
for j=1:numel(ids)
    ix=find(abs(ux-round(x(j),10))<1e-9,1);
    iy=find(abs(uy-round(y(j),10))<1e-9,1);
    xL(j)=xb(ix,1);xR(j)=xb(ix,2);yB(j)=yb(iy,1);yT(j)=yb(iy,2);
end
G=struct('ids',ids,'x',x,'y',y,'xL',xL,'xR',xR,'yB',yB,'yT',yT);
key=newkey;cache=G;
end

function w=footprintWeights(G,xp,yp,hx,hy,a,b,conserve)
% Uniform rectangular rail-to-deck influence kernel integrated over nodal
% tributary areas. Returned weights are dimensionless nodal force fractions.
xl=max(xp-hx,a);xr=min(xp+hx,b);
if xr<=xl,w=[];return;end
ylo=yp-hy;yhi=yp+hy;
dx=max(0,min(G.xR,xr)-max(G.xL,xl));
dy=max(0,min(G.yT,yhi)-max(G.yB,ylo));
area=dx.*dy;
if conserve
    s=sum(area);
    if s<=0,w=[];return;end
    w=area/s;
else
    % Optional transition-footprint mode: only the portion of the nominal
    % footprint overlapping the bridge contributes to the bridge load.
    w=area/(4*hx*hy);
end
w=w(:);
end
