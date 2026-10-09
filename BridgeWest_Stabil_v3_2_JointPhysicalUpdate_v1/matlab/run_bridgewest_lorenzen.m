function fe=run_bridgewest_lorenzen(configFile,outputDir)
% RUN_BRIDGEWEST_LORENZEN Rebuild and solve the provisional Stabil shell model.
% fe=run_bridgewest_lorenzen(); plot_bridgewest_lorenzen(fe);
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1||isempty(configFile),configFile=fullfile(root,'config','model.json');end
if nargin<2||isempty(outputDir),outputDir=fullfile(root,'results','matlab');end
addpath(fullfile(root,'dependencies','Stabil'));c=jsondecode(fileread(configFile));
if ~c.geometry_verified
    warning('BridgeWest:ProvisionalGeometry','Geometry is provisional: see evidence/parameter_provenance.csv. This is not a recovered SOFiSTiK model.');
end
if ~isfolder(outputDir),mkdir(outputDir);end
fe=bw_lorenzen_modes(bw_lorenzen_assemble(bw_lorenzen_mesh(c)));
fe.observation_maps=bw_lorenzen_observation_maps(fe);
fe.modal_table=bw_lorenzen_classify(fe);writetable(fe.modal_table,fullfile(outputDir,'modal_spectrum.csv'));
targets=readtable(fullfile(root,'evidence','thesis_modal_targets.csv'),'TextType','string');
fe.comparison=bw_lorenzen_compare(fe.modal_table,targets);
writetable(fe.comparison,fullfile(outputDir,'thesis_comparison.csv'));
save(fullfile(outputDir,'bridgewest_lorenzen_fe.mat'),'fe','-v7.3');
disp(fe.comparison);fprintf('Nodes: %d; shell elements: %d; free DOFs: %d\n',size(fe.nodes,1),size(fe.shells,1),numel(fe.free));
fprintf('Maximum eigenpair residual: %.3g\n',max(fe.residual));
end
