function idx=bw_find_col(T,aliases,label)
names=string(T.Properties.VariableNames);
key=lower(regexprep(names,'[^A-Za-z0-9]',''));
ak=lower(regexprep(string(aliases),'[^A-Za-z0-9]',''));
idx=find(ismember(key,ak),1);
if isempty(idx)
    % Conservative fallback for common header variants such as
    % "speed [km/h]" or "train_speed_km_h".
    if contains(lower(label),'speed')
        q=contains(key,'speed')|contains(key,'velocity');
    elseif contains(lower(label),'run')
        q=key=="no"|contains(key,'run')|contains(key,'record');
    elseif contains(lower(label),'channel')
        q=contains(key,'channel')&contains(key,'name');
    elseif contains(lower(label),'span')
        q=contains(key,'span')|contains(key,'locationm')|contains(key,'positionm');
    elseif contains(lower(label),'right')
        q=contains(key,'weight')&contains(key,'right');
    elseif contains(lower(label),'left')
        q=contains(key,'weight')&contains(key,'left');
    else
        q=false(size(key));
    end
    idx=find(q,1);
end
if isempty(idx)
    error('BridgeWest:MetadataColumn', ...
        'Could not identify %s in %s. Available variables: %s', ...
        label,inputname(1),strjoin(names,', '));
end
end

