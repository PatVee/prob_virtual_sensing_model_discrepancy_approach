function fe=bw_lorenzen_assemble(mesh)
% Assemble through the supplied Stabil ke_shell4 and ke_beam routines.
% DOF order per node: ux,uy,uz,rx,ry,rz. Units: m,N,kg,s.
c=mesh.c;nodes=mesh.nodes;n=6*size(nodes,1);rho=c.specific_weight_N_m3/c.g_m_s2;
material=[c.E_Pa c.nu rho];
ntrip=576*size(mesh.shells,1)+144*size(mesh.beams,1);
I=zeros(ntrip,1);J=I;Vk=I;Vm=I;pos=0;steelMass=0;
deckArea=zeros(size(nodes,1),1);walkArea=deckArea;foundation=zeros(n,1);
cache=containers.Map('KeyType','char','ValueType','any');
for e=1:size(mesh.shells,1)
    row=mesh.shells(e,:);ids=row(1:4);p=nodes(ids,:);t=row(5);group=row(6);
    relative=p-p(1,:);key=['s' sprintf('%.8g,',[relative(:);t])];
    if isKey(cache,key),km=cache(key);else,[ke,me]=ke_shell4(p,t,material,{});km={ke,me};cache(key)=km;end
    ii=reshape((6*(ids(:)-1)+(1:6))',[],1);add(ii,km{1},km{2});
    area=norm(cross(p(2,:)-p(1,:),p(4,:)-p(1,:)));steelMass=steelMass+rho*t*area;
    if group==4,deckArea(ids)=deckArea(ids)+area/4;end
    if group==3 && (max(p(:,1))<=c.top_patch_length_m+1e-7 || min(p(:,1))>=c.length_m-c.top_patch_length_m-1e-7)
        axis=find('xyz'==c.top_spring_direction,1);dofs=6*(ids-1)+axis;
        foundation(dofs)=foundation(dofs)+c.top_foundation_N_m3*area/4;
    end
end
opt=struct('lumped',false,'rotationalInertia',false);
for e=1:size(mesh.beams,1)
    row=mesh.beams(e,:);ids=row(1:2);p=nodes(ids,:);sec=row(3:6);
    ex=(p(2,:)-p(1,:))/norm(p(2,:)-p(1,:));ref=[0 0 1];
    if abs(ex*ref')>.9,ref=[0 1 0];end
    ey=cross(ref,ex);ey=ey/norm(ey);p3=[p;p(1,:)+ey];
    relative=p-p(1,:);key=['b' sprintf('%.8g,',[relative(:);sec(:)])];
    if isKey(cache,key),km=cache(key);else
        section=[sec(1) Inf Inf sec(2:4)];[ke,me]=ke_beam(p3,section,material,opt);km={ke,me};cache(key)=km;
    end
    ii=reshape((6*(ids(:)-1)+(1:6))',[],1);add(ii,km{1},km{2});
    le=norm(p(2,:)-p(1,:));steelMass=steelMass+rho*sec(1)*le;
    if row(7)==2 && abs(p(2,1)-p(1,1))>1e-7,walkArea(ids)=walkArea(ids)+le*c.walkway_width_m/2;end
end
K=sparse(I(1:pos),J(1:pos),Vk(1:pos),n,n)+spdiags(foundation,0,n,n);
M=sparse(I(1:pos),J(1:pos),Vm(1:pos),n,n);
walkMass=walkArea*c.walkway_surface_mass_kg_m2;
residual=c.mass_per_length_total_kg_m*c.length_m-steelMass-sum(walkMass);
assert(residual>0,'Structural steel/walkway exceed prescribed total mass.');
added=walkMass+residual*deckArea/sum(deckArea);md=zeros(n,1);
for j=1:3,md(j:6:end)=added;end
M=M+spdiags(md,0,n,n);
fe=mesh;fe.K=K;fe.M=M;fe.free=setdiff((1:n)',mesh.fixed);
fe.foundation=foundation;fe.added_mass=added;fe.rho=rho;
fe.steel_mass=steelMass;fe.walkway_mass=sum(walkMass);fe.residual_mass=residual;
    function add(ii,ke,me)
        [ri,ci]=ndgrid(ii,ii);nq=numel(ke);sel=pos+(1:nq);
        I(sel)=ri(:);J(sel)=ci(:);Vk(sel)=ke(:);Vm(sel)=me(:);pos=pos+nq;
    end
end
