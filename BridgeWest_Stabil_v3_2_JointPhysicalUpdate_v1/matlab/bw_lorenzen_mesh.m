function mesh=bw_lorenzen_mesh(c)
% Provisional geometry, with documented source dimensions in config/model.json.
% All node and element connectivity in this MATLAB structure is 1-based.
nodes=zeros(0,3); shells=zeros(0,6); beams=zeros(0,7);
lookup=containers.Map('KeyType','char','ValueType','double');
L=c.length_m; a=c.bearing_x_m(1); b=c.bearing_x_m(2);
ys=c.girder_y_m(:)'; h=c.girder_depth_m; bf=c.flange_width_m; zd=c.deck_z_m;
cross=unique([0 linspace(a,b,round((b-a)/c.crossbeam_spacing_m)+1) L]);
xx=subdivide([0 c.top_patch_length_m L-c.top_patch_length_m L a b cross],c.max_dx_m);
zz=subdivide([0 zd h],c.max_dz_m); yd=subdivide(ys,c.max_dy_m);
for y=ys
    surface(xx,zz,@(x,z)[x y z],c.web_thickness_m,1);
    surface(xx,[y-bf/2 y y+bf/2],@(x,Y)[x Y 0],c.flange_thickness_m,2);
    surface(xx,[y-bf/2 y y+bf/2],@(x,Y)[x Y h],c.flange_thickness_m,3);
end
surface(xx,yd,@(x,y)[x y zd],c.deck_thickness_m,4);
d=c.crossbeam_depth_m; w=c.crossbeam_width_m; tw=c.crossbeam_web_m; tf=c.crossbeam_flange_m;
A=2*w*tf+(d-2*tf)*tw; Iy=(w*d^3-(w-tw)*(d-2*tf)^3)/12;
Iz=2*tf*w^3/12+(d-2*tf)*tw^3/12; J=(2*w*tf^3+(d-2*tf)*tw^3)/3;
section=[A J Iy Iz];
for x=cross
    for iy=1:numel(yd)-1
        beams(end+1,:)=[node([x yd(iy) zd]) node([x yd(iy+1) zd]) section 1]; %#ok<AGROW>
    end
end
if c.include_walkway
    yo=ys(1)-c.walkway_width_m; zr=zd+c.handrail_height_m;
    D=c.handrail_tube_diameter_m; di=D-2*c.handrail_tube_wall_m;
    A=pi*(D^2-di^2)/4; I=pi*(D^4-di^4)/64; rail=[A 2*I I I];
    for x=cross
        beams(end+1,:)=[node([x ys(1) zd]) node([x yo zd]) section 2]; %#ok<AGROW>
        beams(end+1,:)=[node([x yo zd]) node([x yo zr]) rail 3]; %#ok<AGROW>
    end
    for ix=1:numel(cross)-1
        x0=cross(ix);x1=cross(ix+1);
        beams(end+1,:)=[node([x0 yo zd]) node([x1 yo zd]) section 2]; %#ok<AGROW>
        beams(end+1,:)=[node([x0 yo zr]) node([x1 yo zr]) rail 3]; %#ok<AGROW>
    end
end
fixed=[];
for x=[a b]
    for y=ys, fixed(end+1)=6*(node([x y 0])-1)+3; end %#ok<AGROW>
end
n0=node([a ys(1) 0]); n1=node([b ys(1) 0]);
fixed=unique([fixed 6*(n0-1)+[1 2] 6*(n1-1)+2]);
top=zeros(2,numel(xx));
for g=1:2,for ix=1:numel(xx),top(g,ix)=node([xx(ix) ys(g) h]);end,end
mesh=struct('nodes',nodes,'shells',shells,'beams',beams,'fixed',fixed(:),'top',top,'x',xx,'c',c);

    function id=node(p)
        p=round(p*1e8)/1e8;key=sprintf('%.8f,%.8f,%.8f',p);
        if isKey(lookup,key),id=lookup(key);else,id=size(nodes,1)+1;nodes(id,:)=p;lookup(key)=id;end
    end
    function surface(us,vs,fun,t,group)
        grid=zeros(numel(us),numel(vs));
        for i=1:numel(us),for j=1:numel(vs),grid(i,j)=node(fun(us(i),vs(j)));end,end
        for i=1:numel(us)-1
            for j=1:numel(vs)-1
                shells(end+1,:)=[grid(i,j) grid(i+1,j) grid(i+1,j+1) grid(i,j+1) t group]; %#ok<AGROW>
            end
        end
    end
end
function out=subdivide(v,step)
v=unique(round(v*1e9)/1e9);out=[];
for i=1:numel(v)-1
    x=linspace(v(i),v(i+1),ceil((v(i+1)-v(i))/step-1e-8)+1);out=[out x(1:end-1)]; %#ok<AGROW>
end
out=[out v(end)];
end
