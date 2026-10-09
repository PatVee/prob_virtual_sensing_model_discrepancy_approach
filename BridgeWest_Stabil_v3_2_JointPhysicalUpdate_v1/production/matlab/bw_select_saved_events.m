function events=bw_select_saved_events(outputRoot,varargin)
%BW_SELECT_SAVED_EVENTS Select existing 2022/2025 predictions, without raw data.
% Event: one ID, a string/cell array, 'all', '2022', or '2025'.
% Campaign + Runs: all requested run numbers must exist in that campaign.
ip=inputParser;
addParameter(ip,'Event',[]);
addParameter(ip,'Campaign','',@(x)(ischar(x)||isstring(x)||isnumeric(x))&&isscalar(string(x)));
addParameter(ip,'Runs',[],@(x)(isnumeric(x)&&isvector(x))||isempty(x));
parse(ip,varargin{:});o=ip.Results;
assert(isfolder(outputRoot),'BridgeWest:OutputRoot','Result folder does not exist: %s',outputRoot);
d=dir(fullfile(outputRoot,'predictions','*_run*.mat'));avail=strings(0,1);
for k=1:numel(d)
    s=erase(string(d(k).name),'.mat');
    if ~isempty(regexp(char(s),'^202[25]_run\d+$','once')),avail(end+1,1)=s;end %#ok<AGROW>
end
assert(~isempty(avail),'BridgeWest:NoPredictions','No saved 2022/2025 predictions found in %s.',outputRoot);
[yy,rr]=parts(avail);[~,ii]=sortrows([yy rr],[1 2]);avail=avail(ii);yy=yy(ii);rr=rr(ii);
e=string(o.Event);e=e(:);e=e(strlength(e)>0);camp=string(o.Campaign);
if strlength(camp)>0 && ~isempty(e)
    error('BridgeWest:Selection','Use Event OR Campaign/Runs, not both.');
end
if ~isempty(o.Runs)
    assert(strlength(camp)>0,'BridgeWest:Selection','Runs requires Campaign=2022 or 2025.');
    assert(all(isfinite(o.Runs))&&all(o.Runs==fix(o.Runs))&&all(o.Runs>0), ...
        'BridgeWest:Selection','Runs must be positive integer run numbers.');
end
if strlength(camp)>0
    assert(isscalar(camp)&&any(camp==["2022","2025","all"]), ...
        'BridgeWest:Selection','Campaign must be 2022, 2025 or all.');
    if isempty(o.Runs)
        if camp=="all",want=avail;else,want=avail(yy==str2double(camp));end
    else
        assert(camp~="all",'BridgeWest:Selection','Runs requires one campaign, not all.');
        want=strings(numel(o.Runs),1);
        for k=1:numel(o.Runs),want(k)=sprintf('%s_run%02d',char(camp),o.Runs(k));end
    end
elseif isempty(e)
    want="2025_run40";
elseif numel(e)==1 && any(lower(e)==["all","2022","2025"])
    if lower(e)=="all",want=avail;else,want=avail(yy==str2double(e));end
else
    want=strings(numel(e),1);
    for k=1:numel(e)
        q=regexp(char(e(k)),'^(2022|2025)_run(\d+)$','tokens','once');
        assert(~isempty(q),'BridgeWest:Selection','Invalid event ID: %s. Use, for example, 2025_run22.',e(k));
        want(k)=sprintf('%s_run%02d',q{1},str2double(q{2}));
    end
end
want=unique(want(:),'stable');
assert(~isempty(want),'BridgeWest:Selection','No saved predictions match the requested campaign.');
missing=want(~ismember(want,avail));
if ~isempty(missing)
    error('BridgeWest:MissingEvent','No saved prediction for: %s. Available events: %s', ...
        strjoin(cellstr(missing),', '),strjoin(cellstr(avail),', '));
end
events=cellstr(want);
end
function [y,r]=parts(e)
y=zeros(numel(e),1);r=y;
for i=1:numel(e)
    a=regexp(char(e(i)),'^(2022|2025)_run(\d+)$','tokens','once');
    y(i)=str2double(a{1});r(i)=str2double(a{2});
end
end
