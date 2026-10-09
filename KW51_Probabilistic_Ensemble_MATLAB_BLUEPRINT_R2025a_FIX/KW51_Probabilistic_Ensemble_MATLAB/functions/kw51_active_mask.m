function good=kw51_active_mask(E,kin)
t=E.time_s(:);finish=kin(2)+116.316/kin(1)+max(E.rail_offsets_s)+2;
good=t>=kin(2)-1 & t<=finish;
end
