function [B,meta]=bw_stabil_basis_shaker(fe,v,shaker,tobs,c,direction)
%BW_STABIL_BASIS_SHAKER Return strain/displacement unit-bogie response bases.
% B{1}: T x 7 x 8 microstrain/N
% B{2}: T x 41 x 8 mm/N, girder BA
% B{3}: T x 41 x 8 mm/N, girder BB
if nargin<6 || isempty(direction),direction=1;end
if ~isfield(fe,'observation_maps') || isempty(fe.observation_maps)
    fe.observation_maps=bw_lorenzen_observation_maps(fe);
end
[ts,H,state]=bw_stabil_sim_shaker(fe,v,shaker,c,direction);sel=state.selected_modes;
Sm=fe.observation_maps.strain_modal(:,sel);Dm=fe.observation_maps.displacement_modal(:,sel);
if nargin<4 || isempty(tobs),tobs=ts;end
tobs=tobs(:);B={zeros(numel(tobs),7,8),zeros(numel(tobs),41,8),zeros(numel(tobs),41,8)};
for ib=1:8
    q=squeeze(H(:,:,ib));
    ys=(q*Sm.')*1e6; yd=(q*Dm.')*1e3;
    B{1}(:,:,ib)=interp1(ts,ys,tobs,'linear',0);
    B{2}(:,:,ib)=interp1(ts,yd(:,1:41),tobs,'linear',0);
    B{3}(:,:,ib)=interp1(ts,yd(:,42:82),tobs,'linear',0);
end
meta=struct('simulation_time_s',ts,'observation_time_s',tobs,'modal_state',state, ...
    'strain_channels',{fe.observation_maps.strain_channels},'strain_unit','microstrain per N', ...
    'displacement_unit','mm per N','load_basis',loadBasisDescription(c));
end

function s=loadBasisDescription(c)
if isfield(c,'load_transfer_model') && contains(lower(string(c.load_transfer_model)),'two_rail')
    s='8 bogies; 2 axles/bogie; 2 wheels/axle; two-rail distributed deck footprint';
else
    s='8 bogies; two axles per bogie; track-centre resultant';
end
end
