function z=bw_rail_obj(x,t,a,ax,on,v0),if x(2)<.85*v0||x(2)>1.15*v0||abs(x(1)-on)>.35,z=1e12;else,z=-sum(interp1(t,a,x(1)+ax/x(2),'linear',0));end,end
