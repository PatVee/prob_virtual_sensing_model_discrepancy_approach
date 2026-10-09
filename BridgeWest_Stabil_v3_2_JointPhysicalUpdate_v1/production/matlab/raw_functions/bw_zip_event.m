function p=bw_zip_event(A,year,run,isPLS)
if nargin<4 || ~isPLS, prefix='LTW'; else, prefix='PLSW'; end
% [.] avoids sprintf interpreting a backslash escape before regex sees it.
pat=sprintf('^%s_ICETD_%02d_DT_%d_.*[.]h5$',prefix,run,year);
hit=find(~cellfun(@isempty,regexpi({A.members.name},pat)));
assert(numel(hit)==1,'BridgeWest:MissingEvent','Expected one %s record for %d run %d; found %d.',prefix,year,run,numel(hit));
p=bw_zip_get(A,A.members(hit).name);
end
