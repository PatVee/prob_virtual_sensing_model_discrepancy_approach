function bw_verify_bridgewest_sensor_coordinates(sensorFile,names)
%BW_VERIFY_BRIDGEWEST_SENSOR_COORDINATES Guard the Bridge West channel mapping.
S=bw_read_csv(sensorFile);
ich=bw_find_col(S,["channelname","channel","name"],"channel name");
ix=bw_find_col(S,["spanm","span","locationm","xpositionm"],"span position");
sch=bw_strcol(S,ich);sx=bw_numcol(S,ix);
for i=1:numel(names)
    ch=char(string(names{i}));
    switch ch
        case {'str_x_BA_uf_031','str_x_BB_uf_031'}
            xexp=5.62;
        case {'str_x_BA_uf_053','str_x_BB_uf_053'}
            xexp=9.72;
        case 'str_x_BA_uf_058'
            xexp=10.74;
        case 'str_x_BA_uf_070'
            xexp=12.80;
        case 'str_x_BB_uf_081'
            xexp=14.85;
        otherwise
            error('BridgeWest:UnexpectedChannel','No verified Bridge West coordinate is registered for %s.',ch);
    end
    j=find(sch==string(ch),1);
    assert(~isempty(j),'Bridge West sensor metadata does not contain %s.',ch);
    assert(abs(sx(j)-xexp)<=0.05,'Unexpected Bridge West coordinate for %s: %.4f m.',ch,sx(j));
end
end
