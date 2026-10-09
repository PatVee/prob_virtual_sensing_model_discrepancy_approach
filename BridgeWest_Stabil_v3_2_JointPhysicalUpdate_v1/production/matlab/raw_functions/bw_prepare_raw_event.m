function z=bw_prepare_raw_event(A,eid,c,outputRoot,names)
%BW_PREPARE_RAW_EVENT Prepare ONLY the explicitly permitted inference strains.
% In production v3, NAMES is A031/B053 only. Protected target channels are
% never read or preprocessed by this function.
if nargin<5||isempty(names)
    error('BridgeWest:ExplicitChannelsRequired', ...
        'Inference preprocessing requires an explicit permitted-channel list.');
end
names=cellstr(string(names(:)).');
tok=regexp(eid,'^(\d{4})_run(\d+)$','tokens','once');assert(~isempty(tok),'Invalid event ID');
year=str2double(tok{1});run=str2double(tok{2});
if year==2022,suffix='21';else,suffix='24';end
sensor=bw_zip_get(A,['sensors_bridge_LT_W' suffix '.csv']);
protocol=bw_zip_get(A,['protocol_WTD_' sprintf('%02d',mod(year,100)) '.csv']);
path=bw_zip_event(A,year,run,false);
sp=bw_protocol_speeds(protocol);assert(run<=numel(sp)&&isfinite(sp(run))&&sp(run)>0,'Invalid protocol speed.');
[v,entry]=bw_kinematics(path,sensor,sp(run),c);
bw_verify_bridgewest_sensor_coordinates(sensor,names);
m=bw_prepstrain(path,names,v,entry,c);
z.time=m.time;z.speed=v;z.entry=entry;z.channels={};z.flags={};z.Y=[];z.source=path;
for i=1:numel(names)
    ch=names{i};nm=matlab.lang.makeValidName(ch);
    if ~isfield(m.Y,nm),continue;end
    z.channels{end+1}=ch;z.flags{end+1}=char(m.flags.(nm));z.Y(:,end+1)=m.Y.(nm); %#ok<AGROW>
end
assert(numel(z.channels)==numel(names),'Not all explicitly permitted inference channels were prepared.');
processed=fullfile(outputRoot,'processed_references');if ~exist(processed,'dir'),mkdir(processed);end
save(fullfile(processed,[eid '.mat']),'-struct','z','-v7');
pr=struct('event',eid,'source',path,'source_sha256',bw_sha256(path), ...
    'speed_m_s',v,'entry_s',entry,'timing','Upstream rail strain only', ...
    'sensor_metadata',sensor,'protocol',protocol,'preprocessed_channels',{names}, ...
    'protected_target_histories_loaded',false);
bw_write_json(fullfile(processed,[eid '_source.json']),pr);
end
