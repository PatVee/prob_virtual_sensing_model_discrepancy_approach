function maps=bw_lorenzen_observation_maps(fe)
% Longitudinal surface strain from Stabil shell stresses, not a fitted depth.
% x_model=18.4-x_dataset. Gauge transverse offsets are nominal (outer flange).
% Confirm gauge side and plate thicknesses before using for field prediction.
names={'A031','B031','A053_lower','B053','A058','A070','B081'};
xd=[5.62 5.62 9.72 9.72 10.74 12.8 14.85];girder=[1 2 1 2 1 1 2];
nd=size(fe.K,1);Bstrain=sparse(7,nd);E=fe.c.E_Pa;nu=fe.c.nu;
material=[E nu fe.rho];
for k=1:7
    x=fe.c.length_m-xd(k);g=girder(k);signSide=2*g-3;
    y=fe.c.girder_y_m(g)+signSide*fe.c.flange_width_m/4;
    group=3;offset=48;if k==3,group=2;offset=0;end
    candidates=[];
    for e=find(fe.shells(:,6)==group)'
        p=fe.nodes(fe.shells(e,1:4),:);
        if x>=min(p(:,1))-1e-8 && x<=max(p(:,1))+1e-8 && y>=min(p(:,2))-1e-8 && y<=max(p(:,2))+1e-8
            candidates(end+1)=e; %#ok<AGROW>
        end
    end
    assert(~isempty(candidates),'Gauge falls outside flange mesh.');
    for e=candidates
        ids=fe.shells(e,1:4);p=fe.nodes(ids,:);t=fe.shells(e,5);
        xi=2*(x-min(p(:,1)))/(max(p(:,1))-min(p(:,1)))-1;
        eta=2*(y-min(p(:,2)))/(max(p(:,2))-min(p(:,2)))-1;
        N=.25*[(1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];
        stress=se_shell4(p,t,material,eye(24),{},'cart');
        strain=zeros(4,24);
        for j=1:4
            ii=offset+6*(j-1);strain(j,:)=(stress(ii+1,:)-nu*(stress(ii+2,:)+stress(ii+3,:)))/E;
        end
        dof=reshape((6*(ids(:)-1)+(1:6))',[],1);
        Bstrain(k,dof)=Bstrain(k,dof)+N*strain/numel(candidates);
    end
end
xs=linspace(fe.c.bearing_x_m(1),fe.c.bearing_x_m(2),41);Bdisp=sparse(82,nd);
for g=1:2
    for j=1:numel(xs)
        k=find(fe.x<=xs(j)+1e-9,1,'last');k=min(k,numel(fe.x)-1);
        a=(xs(j)-fe.x(k))/(fe.x(k+1)-fe.x(k));dof=6*(fe.top(g,k:k+1)-1)+3;
        Bdisp((g-1)*41+j,dof)=[1-a a];
    end
end
maps=struct('strain',Bstrain,'displacement',Bdisp,'strain_modal',Bstrain*fe.P, ...
    'displacement_modal',Bdisp*fe.P,'strain_channels',{names},'x_model_m',xs, ...
    'x_dataset_m',fe.c.length_m-xs,'strain_units','dimensionless','displacement_units','m', ...
    'provisional_gauge_transverse_offset',true);
end
