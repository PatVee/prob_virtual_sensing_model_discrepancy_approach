function [data,names,rates,flags,fsr]=bw_readstd_selected(file,requested)
%BW_READSTD_SELECTED Read only selected HDF5 signal columns.
% Metadata are read in full because they contain channel labels/rates/flags;
% signal histories are read only for REQUESTED channels.
md=string(h5read(file,'/metadata').');allNames=md(1,:);allRates=str2double(md(3,:));allFlags=md(4,:);
fsr=max(allRates);
requested=string(requested(:)).';idx=zeros(1,numel(requested));
for k=1:numel(requested)
    j=find(allNames==requested(k));assert(numel(j)==1,'Requested HDF5 channel %s missing or duplicated.',requested(k));idx(k)=j;
end
names=allNames(idx);rates=allRates(idx);flags=allFlags(idx);
info=h5info(file,'/data');sz=double(info.Dataspace.Size);nch=numel(allNames);
assert(numel(sz)==2,'Expected a 2-D /data dataset.');
if sz(1)==nch
    ns=sz(2);data=nan(ns,numel(idx));
    for k=1:numel(idx)
        x=h5read(file,'/data',[idx(k) 1],[1 ns]);data(:,k)=double(x(:));
    end
elseif sz(2)==nch
    ns=sz(1);data=nan(ns,numel(idx));
    for k=1:numel(idx)
        x=h5read(file,'/data',[1 idx(k)],[ns 1]);data(:,k)=double(x(:));
    end
else
    error('BridgeWest:HDF5Layout','Neither dimension of /data matches the metadata channel count.');
end
end
