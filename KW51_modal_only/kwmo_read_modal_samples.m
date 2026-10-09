function S=kwmo_read_modal_samples(csvFile,bank)
% Portable reader for the generated numeric modal CSV; timestamp IDs required.
fid=fopen(csvFile,'r');assert(fid>=0,'Cannot read modal CSV.');clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
header=strsplit(strrep(strtrim(fgetl(fid)),'"',''),',');
required={'sample_id','fe_mode','frequency_Hz','damping_ratio'};
[ok,ix]=ismember(required,header);assert(all(ok),'Modal CSV is missing required columns.');
ids={};mode=[];freq=[];damp=[];
while true
 line=fgetl(fid);if ~ischar(line),break;end;if isempty(strtrim(line)),continue;end
 cells=strsplit(strrep(line,'"',''),',','CollapseDelimiters',false);
 assert(numel(cells)>=max(ix),'Malformed modal CSV row.');
 ids{end+1,1}=strtrim(cells{ix(1)}); %#ok<AGROW>
 mode(end+1,1)=str2double(cells{ix(2)});freq(end+1,1)=str2double(cells{ix(3)});damp(end+1,1)=str2double(cells{ix(4)}); %#ok<AGROW>
end
samples=unique(ids,'stable');modes=unique(mode,'sorted');
assert(numel(samples)>=24,'At least 24 complete measured joint samples required.');
assert(all(isfinite(modes)&modes>=1&modes<=numel(bank.frequency_Hz)&modes==floor(modes)),'Invalid FE mode indices.');
f=zeros(numel(samples),numel(modes));z=f;
for i=1:numel(samples)
 for j=1:numel(modes)
  k=strcmp(ids,samples{i})&mode==modes(j);
  assert(sum(k)==1,'Each sample must contain every selected FE mode exactly once.');
  f(i,j)=freq(k);z(i,j)=damp(k);
 end
end
assert(all(isfinite(f(:))&f(:)>0),'Invalid frequencies.');
assert(all(isfinite(z(:))&z(:)>0&z(:)<.1),'Expected screened damping ratios, not percent.');
sdn=datenum(samples,'yyyymmdd_HHMMSS');
S=struct('sample_ids',{samples},'sample_sdn',sdn,'fe_modes',modes(:), ...
 'frequency_Hz',f,'damping_ratio',z,'frequency_ratio',f./median(f,1), ...
 'median_damping',median(z,1),'source_sha256',kw51_sha256(csvFile), ...
 'policy','Joint empirical modal resampling; frequency ratios centred on the unscaled FE; no target calibration.');
end
