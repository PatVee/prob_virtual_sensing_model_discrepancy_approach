function F=bw_stabil_force_map(fe,t,v,c,direction)
%BW_STABIL_FORCE_MAP Dispatch the configured moving-load transfer operator.
if nargin<5 || isempty(direction),direction=1;end
model='centerline_point';
if isfield(c,'load_transfer_model') && ~isempty(c.load_transfer_model)
    model=lower(strtrim(char(c.load_transfer_model)));
end
switch model
    case {'centerline_point','centreline_point','legacy_centerline'}
        F=bw_stabil_force_centerline(fe,t,v,c,direction);
    case {'two_rail_distributed','two_rail_ballast','two_rail'}
        F=bw_stabil_force_two_rail(fe,t,v,c,direction);
    otherwise
        error('BridgeWest:UnknownLoadTransfer','Unknown load_transfer_model: %s',model);
end
end
