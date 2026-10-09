function z=bw_prepare_validation_targets(A,eid,c,outputRoot,names,v,entry,expectedTime)
%BW_PREPARE_VALIDATION_TARGETS Read protected targets only AFTER prediction.
% Kinematics are inherited from the reference-only prediction stage; target
% measurements cannot change speed, entry time, loads, discrepancy or FE state.
names=cellstr(string(names(:)).');
tok=regexp(eid,'^(\d{4})_run(\d+)$','tokens','once');assert(~isempty(tok),'Invalid event ID');
year=str2double(tok{1});run=str2double(tok{2});
if year==2022,suffix='21';else,suffix='24';end
sensor=bw_zip_get(A,['sensors_bridge_LT_W' suffix '.csv']);
path=bw_zip_event(A,year,run,false);
bw_verify_bridgewest_sensor_coordinates(sensor,names);
m=bw_prepstrain(path,names,v,entry,c);
assert(numel(m.time)==numel(expectedTime),'Validation target time grid differs from frozen prediction grid.');
assert(max(abs(double(m.time(:))-double(expectedTime(:))))<1e-10,'Validation target time grid is not identical to prediction grid.');
z.time=m.time;z.speed=v;z.entry=entry;z.channels={};z.flags={};z.Y=[];z.source=path;
for i=1:numel(names)
    ch=names{i};nm=matlab.lang.makeValidName(ch);
    if ~isfield(m.Y,nm)
        z.channels{end+1}=ch;z.flags{end+1}='missing';z.Y(:,end+1)=NaN(size(m.time)); %#ok<AGROW>
    else
        z.channels{end+1}=ch;z.flags{end+1}=char(m.flags.(nm));z.Y(:,end+1)=m.Y.(nm); %#ok<AGROW>
    end
end
out=fullfile(outputRoot,'validation_targets');if ~exist(out,'dir'),mkdir(out);end
save(fullfile(out,[eid '.mat']),'-struct','z','-v7');
meta=struct('event',eid,'source',path,'source_sha256',bw_sha256(path), ...
    'speed_m_s',v,'entry_s',entry,'validation_channels',{names}, ...
    'stage','assessment only; created after prediction');
bw_write_json(fullfile(out,[eid '_source.json']),meta);
end
