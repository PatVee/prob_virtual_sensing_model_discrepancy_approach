function [data,names,rates,flags]=bw_readstd(file),data=h5read(file,'/data').';md=string(h5read(file,'/metadata').');names=md(1,:);rates=str2double(md(3,:));flags=md(4,:);end
