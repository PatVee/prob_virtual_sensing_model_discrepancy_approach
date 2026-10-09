function inputs=kw51_prepare(catalog,inputFolder,events,maxEvents)
if ~exist(inputFolder,'dir'),mkdir(inputFolder);end
refs={'sgBD1011A','sgBD1415A','sgBD2728A'};rails={'sgRAAB','sgRAAL','sgRABB','sgRABL'};
inputs=struct('event_id',{},'input_file',{},'raw_file',{});audit={};
if ischar(events),events={events};elseif isstring(events),events=cellstr(events);end
for k=1:numel(catalog)
 C=catalog(k);if ~any(strcmp(events,'all'))&&~any(strcmp(events,C.event_id)),continue;end
 raw=load(C.path,'predat_sg');reason='';ratio=NaN;count=0;
 if ~isfield(raw,'predat_sg'),audit(end+1,:)={C.event_id,0,ratio,count,'Missing predat_sg'};continue;end %#ok<AGROW>
 s=raw.predat_sg;names=cellstr(s.labels);[ok,ix]=ismember([refs rails],names);
 if ~all(ok),audit(end+1,:)={C.event_id,0,ratio,count,'Missing reference/rail channel'};continue;end %#ok<AGROW>
 fs=double(s.fs);a=double(s.tdata(:,ix))*1e6;clear raw
 if any(any(~isfinite(a(:,1:3)))),audit(end+1,:)={C.event_id,0,ratio,count,'Nonfinite reference data'};continue;end %#ok<AGROW>
 for j=4:7,if any(~isfinite(a(:,j))),a(:,j)=0;end;end
 a=a-median(a(1:min(floor(fs),size(a,1)),:),1);tt=(0:size(a,1)-1).'/fs;
 y=kw51_filter(a(:,1:3),10,fs);rail=kw51_filter(a(:,4:7),40,fs);
 tr=(0:ceil(tt(end)/.005)-1).'* .005;rail=interp1(tt,rail,tr);
 scale=median(abs(rail-median(rail,1)),1);energy=sum((rail./max(scale,1e-9)).^2,1);
 te=[sum(energy(1:2)),sum(energy(3:4))];[~,which]=max(te);ratio=max(te)/max(min(te),1e-9);
 cols=(which-1)*2+(1:2);[~,j]=max(sum(rail(:,cols).^2,1));j=cols(j);
 env=conv(abs(rail(:,j)),ones(3,1)/3,'same');med=median(env);sd=1.4826*median(abs(env-med));
 peaks=kw51_peaks(env,max(med+4*sd,.08*max(env)),12);count=numel(peaks);
 accepted=ratio>=1.5 && count>=6 && count<=40 && all(max(y)-min(y)>1);
 if ~accepted,reason='Reference/rail quality gate';end
 audit(end+1,:)={C.event_id,accepted,ratio,count,reason}; %#ok<AGROW>
 if ~accepted,continue;end
 t=(0:ceil(min(tt(end),35)/.02)-1).'* .02;
 E=struct();E.time_s=t;E.y_ref=interp1(tt,y,t);E.reference_names=refs;
 E.rail_offsets_s=tr(peaks)-tr(peaks(1));E.rail_first_time_s=tr(peaks(1));
 E.rail_relative_load=env(peaks)/median(env(peaks));tracks='AB';E.track=tracks(which);
 E.event_id=C.event_id;E.source_member=C.source_member;E.strain_start_sdn=double(s.sdn(1));E.raw_archive=C.raw_archive;
 kw51_check_input(E);path=fullfile(inputFolder,[E.event_id '_input.mat']);save(path,'-struct','E','-v7');
 inputs(end+1)=struct('event_id',E.event_id,'input_file',path,'raw_file',C.path); %#ok<AGROW>
 fprintf('Prepared %s: track %s, %d rail peaks\n',E.event_id,E.track,count);
 if numel(inputs)>=maxEvents,break;end
end
kw51_write_csv(fullfile(inputFolder,'selection_audit.csv'),{'event','accepted','rail_ratio','rail_peaks','reason'},audit);
assert(~isempty(inputs),'No passages passed the permitted-reference and rail checks.');
end
