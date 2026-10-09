function out=test_two_rail_force_map(packageRoot)
%TEST_TWO_RAIL_FORCE_MAP Mechanics-only verification of the new force map.
if nargin<1 || isempty(packageRoot),packageRoot=fileparts(fileparts(mfilename('fullpath')));end
addpath(fullfile(packageRoot,'matlab'));
addpath(fullfile(packageRoot,'dependencies','Stabil'));
addpath(fullfile(packageRoot,'shaker_informed','matlab'));
c=jsondecode(fileread(fullfile(packageRoot,'shaker_informed','discrepancy_config.json')));
fe=run_bridgewest_lorenzen(fullfile(packageRoot,'config','model.json'),fullfile(tempdir,'bw_two_rail_test_fe'));

% Choose bogie 1 with both axles well inside the bearing span.
v=20; a=fe.c.bearing_x_m(1); xlead=a+5.0; t=(xlead-a)/v;
F=bw_stabil_force_two_rail(fe,t,v,c,1);
vertical=3:6:size(F,1);
q=F(vertical,1);
ids=(vertical(:)-3)/6+1;
y=fe.nodes(ids,2);

resultant=-sum(q);
ycent=sum((-q).*y)/sum(-q);
ytrack=mean(fe.c.girder_y_m);
left=-sum(q(y<ytrack));right=-sum(q(y>ytrack));

assert(abs(resultant-1)<1e-10,'Unit bogie resultant not conserved.');
assert(abs(ycent-ytrack)<1e-10,'Two-rail map is not transversely symmetric.');
assert(abs(left-right)<1e-10,'Left/right rail load is not balanced.');

out=struct('unit_bogie_resultant_N',resultant,'transverse_centroid_m',ycent, ...
    'track_centroid_m',ytrack,'left_half_resultant_N',left,'right_half_resultant_N',right, ...
    'status','PASS');
fprintf('Two-rail force-map test PASS: resultant %.12g N, y-centroid %.6f m.\n',resultant,ycent);
end
