function [v,entry]=bw_kinematics(file,sensorFile,vguess,cfg)
%BW_KINEMATICS Infer speed/entry from upstream rail strains only.
% Structural bridge target histories are not read by this function.
S=bw_read_csv(sensorFile);
iCh=bw_find_col(S,["channelname","channel","name"],"channel name");
iX=bw_find_col(S,["spanm","span","locationm","xpositionm"],"span position");
sch=bw_strcol(S,iCh);sx=bw_numcol(S,iX);
% Read HDF5 metadata only to identify available upstream rail timing channels.
md=string(h5read(file,'/metadata').');available=md(1,:);
railNames=strings(0,1);railLoc=[];
for i=1:height(S)
    ch=sch(i);loc=sx(i);
    if startsWith(ch,'str_T')&&isfinite(loc)&&loc>17.4&&any(available==ch)
        railNames(end+1,1)=ch;railLoc(end+1,1)=loc; %#ok<AGROW>
    end
end
assert(~isempty(railNames),'No upstream rail timing channels are available.');
[data,names,rates,~,fsr]=bw_readstd_selected(file,railNames);
ax=double(cfg.axle_offsets_m(:));cand=[];
for i=1:numel(railNames)
    ch=railNames(i);loc=railLoc(i);j=find(names==ch,1);col=data(:,j);ii=find(isfinite(col));
    if numel(ii)<1000,continue;end
    t=double(ii-1)/fsr;y=double(col(ii));fs=rates(j);
    y=y-median(y(1:min(numel(y),round(fs))));y=bw_lowpass(y,fs,min(100,.4*fs));
    ab=abs(y);thr=max(10,.2*bw_percentile(ab,99.5));qq=find(t>2&ab>thr,1);if isempty(qq),continue;end
    on=t(qq);v0=vguess/3.6;tg=linspace(on-.15,on+.25,401);sc=zeros(size(tg));
    for k=1:numel(tg),sc(k)=sum(interp1(t,ab,tg(k)+ax/v0,'linear',0));end
    [~,im]=max(sc);z0=[tg(im) v0];obj=@(z)bw_rail_obj(z,t,ab,ax,on,v0);
    z=fminsearch(obj,z0,optimset('Display','off','MaxIter',200,'TolX',1e-7));
    cand=[cand;z(2),z(1)+(loc-17.4)/z(2)]; %#ok<AGROW>
end
assert(~isempty(cand),'No valid upstream rail timing solution.');
mv=median(cand(:,1));cand=cand(abs(cand(:,1)-mv)/mv<.03,:);me=median(cand(:,2));
mdv=median(abs(cand(:,2)-me))+eps;keep=abs(cand(:,2)-me)<max(.05,3*1.4826*mdv);if any(keep),cand=cand(keep,:);end
v=median(cand(:,1));entry=median(cand(:,2));
end
