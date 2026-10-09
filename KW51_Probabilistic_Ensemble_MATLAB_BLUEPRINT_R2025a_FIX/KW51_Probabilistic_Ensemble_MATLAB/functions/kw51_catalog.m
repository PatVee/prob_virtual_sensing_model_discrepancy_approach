function catalog=kw51_catalog(dataSource,cacheRoot)
if ischar(dataSource),dataSource={dataSource};elseif isstring(dataSource),dataSource=cellstr(dataSource);end
assert(iscell(dataSource)&&~isempty(dataSource),'Provide a ZIP path, folder, or cell array of ZIP paths.');
if ~exist(cacheRoot,'dir'),mkdir(cacheRoot);end
raw={};archives={};
for j=1:numel(dataSource)
 p=char(dataSource{j});assert(exist(p,'file')||exist(p,'dir'),'Data source does not exist: %s',p);
 if exist(p,'dir')
  archives=[archives;kw51_files(p,'(?i)\.zip$')]; %#ok<AGROW>
  raw=[raw;kw51_files(p,'^traindata_[0-9]{8}_[0-9]{6}\.mat$')]; %#ok<AGROW>
 else
  [~,~,ext]=fileparts(p);
  if strcmpi(ext,'.zip'),archives{end+1,1}=p;elseif strcmpi(ext,'.mat'),raw{end+1,1}=p;else,error('Unsupported data file: %s',p);end
 end
end
catalog=struct('event_id',{},'path',{},'source_member',{},'raw_archive',{});
for j=1:numel(raw)
 [~,event]=fileparts(raw{j});catalog(end+1)=struct('event_id',event,'path',raw{j},'source_member',raw{j},'raw_archive',''); %#ok<AGROW>
end
archives=unique(archives);
for j=1:numel(archives)
 z=archives{j};fprintf('Indexing ZIP: %s\n',z);h=kw51_sha256(z);
 target=fullfile(cacheRoot,['zip_' h(1:16)]);marker=fullfile(target,'extraction_complete.mat');
 if ~exist(marker,'file')
  if ~exist(target,'dir'),mkdir(target);end
  unzip(z,target);archive_sha256=h;save(marker,'archive_sha256','-v7');
 else
  m=load(marker);assert(strcmp(m.archive_sha256,h),'Archive cache hash mismatch.');
 end
 members=kw51_files(target,'^traindata_[0-9]{8}_[0-9]{6}\.mat$');
 for k=1:numel(members)
  [~,event]=fileparts(members{k});rel=members{k}(numel(target)+2:end);
  catalog(end+1)=struct('event_id',event,'path',members{k},'source_member',rel,'raw_archive',z); %#ok<AGROW>
 end
end
assert(~isempty(catalog),'No raw traindata_YYYYMMDD_HHMMSS.mat records found.');
[~,ix]=sort({catalog.event_id});catalog=catalog(ix);keep=true(size(catalog));
for k=2:numel(catalog)
 if strcmp(catalog(k).event_id,catalog(k-1).event_id)
  assert(strcmp(kw51_sha256(catalog(k).path),kw51_sha256(catalog(k-1).path)), ...
   'Conflicting files have the same event ID: %s',catalog(k).event_id);keep(k)=false;
 end
end
catalog=catalog(keep);
end
