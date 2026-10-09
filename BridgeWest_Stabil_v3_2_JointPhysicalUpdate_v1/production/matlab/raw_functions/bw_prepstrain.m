function m=bw_prepstrain(file,chs,v,entry,cfg)
[data,names,rates,flags,fsr]=bw_readstd_selected(file,chs); pre=cfg.pre_entry_s; post=cfg.post_exit_s;
te=entry+(cfg.free_span_m+max(cfg.axle_offsets_m))/v+post;
tg=(entry-pre:1/cfg.fs_output_Hz:te).';
m.time=tg-entry; m.speed=v; m.Y=struct(); m.flags=struct();
for ch=string(chs(:)).'
 j=find(names==ch,1); if isempty(j), continue; end
 col=data(:,j); ii=find(isfinite(col)); nm=matlab.lang.makeValidName(ch);
 % An unusable target channel must not stop reference-based prediction.
 % Fourth-order filtfilt requires more than 12 input samples.
 if numel(ii)<=12 || ~isfinite(rates(j)) || rates(j)<=2*cfg.strain_lowpass_Hz
  m.Y.(nm)=NaN(size(tg)); m.flags.(nm)='insufficient_data'; continue;
 end
 t=double(ii-1)/fsr; y=double(col(ii)); q=t>=entry-pre & t<entry;
 if any(q), base=median(y(q));
 else, base=median(y(1:min(numel(y),round(rates(j))))); end
 y=bw_lowpass(y-base,rates(j),cfg.strain_lowpass_Hz);
 m.Y.(nm)=interp1(t,y,tg,'linear','extrap'); m.flags.(nm)=flags(j);
end
end
