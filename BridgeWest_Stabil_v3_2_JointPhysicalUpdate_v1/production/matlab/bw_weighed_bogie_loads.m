function W=bw_weighed_bogie_loads(campaign,weightDir)
%BW_WEIGHED_BOGIE_LOADS Read original campaign static weighing, not FE estimates.
% Returns eight bogie loads in protocol order. Static weighing is NOT a
% measurement of time-varying wheel/rail contact forces during a passage.
if nargin<2||isempty(weightDir)
    here=fileparts(mfilename('fullpath'));pkg=fileparts(fileparts(here));
    weightDir=fullfile(pkg,'data','weighing');
end
year=str2double(string(campaign));
assert(isscalar(year)&&ismember(year,[2022 2025]),'BridgeWest:WeighingCampaign','Campaign must be 2022 or 2025.');
source=fullfile(weightDir,sprintf('weights_ICE_TD_%d.csv',year));
assert(isfile(source),'BridgeWest:WeighingFile','Original weighing CSV is missing: %s',source);
T=readtable(source);required={'no','location_m','weight_right_kg','weight_left_kg'};
assert(all(ismember(required,T.Properties.VariableNames)),'BridgeWest:WeighingSchema','Unexpected weighing columns.');
T=sortrows(T,'no');
assert(height(T)==16&&isequal(double(T.no(:)),(1:16).'),'BridgeWest:WeighingSchema','Expected 16 ordered axle stations.');
expected=[0 2.6 19 21.6 25.9 28.5 44.9 47.5 51.8 54.4 70.8 73.4 77.7 80.3 96.7 99.3].';
assert(max(abs(double(T.location_m(:))-expected))<1e-8,'BridgeWest:WeighingGeometry','Weighing axle offsets differ from the model train.');
a=double(T.weight_right_kg)+double(T.weight_left_kg);
assert(all(isfinite(a))&&all(a>0),'BridgeWest:WeighingMass','Invalid axle masses.');
g=9.80665;p=sum(reshape(a,2,8),1)*g/1000;
W=struct('campaign',year,'source_file',source,'gravity_m_s2',g, ...
    'axle_offsets_m',expected,'axle_mass_kg',a,'bogie_load_kN',p, ...
    'total_mass_kg',sum(a),'total_load_kN',sum(p), ...
    'kind','independent static campaign weighing; not measured dynamic contact forces');
end
