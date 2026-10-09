function truth=bw_pls_raw_snapshot(file)
%BW_PLS_RAW_SNAPSHOT Assessment-only PLS own-peak displacement profile.
%
% This function is never called during train-load/discrepancy inference. It
% reads the PLS file only after the passage prediction has already been saved.
% For each girder, the measured spatial profile is taken at the time of that
% girder's measured displacement peak near bridge midspan.
md=string(h5read(file,'/metadata').');
axesNames={'BA','BB'};
rows=struct('axis',{},'x_m',{},'measured_mm',{},'method',{});
info=h5info(file); datasets=string({info.Datasets.Name});
for ig=1:2
    g=axesNames{ig}; dn=['disp_data_' g]; tn=['time_data_' g];
    if ~any(datasets==string(dn)), continue; end
    D=double(h5read(file,['/' dn]).');
    TT=double(h5read(file,['/' tn]).');
    ids=find(md(2,:)==string(dn));
    x=str2double(md(8,ids)); flags=md(7,ids);
    t=mean(TT,2,'omitnan');
    dt=median(diff(t),'omitnan'); assert(isfinite(dt)&&dt>0,'Invalid PLS time base in %s.',file);
    fs=1/dt;
    good=find(flags=="ok" & isfinite(x) & x>=1 & x<=17.4);
    cols=[]; xx=[];
    % D columns are local to this disp_data_* dataset; metadata columns are global.
    for kk=good(:).'
        y=D(:,kk); valid=isfinite(y); if nnz(valid)<20, continue; end
        ii=find(valid); allidx=(1:numel(y)).';
        y=interp1(ii,y(valid),allidx,'linear','extrap');
        y(allidx<ii(1))=D(ii(1),kk); y(allidx>ii(end))=D(ii(end),kk);
        nbase=max(5,min(numel(y),floor(fs)));
        y=y-mean(y(1:nbase),'omitnan');
        fc=min(10,0.4*fs);
        if fc>0 && fc<fs/2
            [b,a]=butter(4,fc/(fs/2),'low'); y=filtfilt(b,a,y);
        end
        cols(:,end+1)=y; xx(end+1)=x(kk); %#ok<AGROW>
    end
    if isempty(cols), continue; end
    mid=mean([1.0 17.4]); [~,jmid]=min(abs(xx-mid));
    [~,it]=max(abs(cols(:,jmid)));
    peak=cols(it,:);
    for j=1:numel(xx)
        rows(end+1)=struct('axis',g,'x_m',xx(j),'measured_mm',peak(j), ... %#ok<AGROW>
            'method','Measured PLS own-midspan-peak profile');
    end
end
assert(~isempty(rows),'No valid PLS displacement channels in %s.',file);
truth=struct2table(rows);
end
