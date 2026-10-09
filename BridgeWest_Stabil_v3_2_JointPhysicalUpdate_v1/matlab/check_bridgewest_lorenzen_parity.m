function report=check_bridgewest_lorenzen_parity()
% Run locally in MATLAB to compare assembled Stabil matrices with Python.
% This check was supplied, not executed, in the hosted Python environment.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'dependencies','Stabil'));
py=load(fullfile(root,'results','python_model.mat'));c=jsondecode(py.config_json);
fe=bw_lorenzen_assemble(bw_lorenzen_mesh(c));
assert(isequal(size(fe.nodes),size(py.nodes)) && max(abs(fe.nodes-py.nodes),[],'all')<1e-7,'Mesh differs.');
assert(isequal(fe.shells(:,1:4),py.shells(:,1:4)+1),'Connectivity differs.');
report.stiffness_relative_error=norm(fe.K-py.K,'fro')/norm(py.K,'fro');
report.mass_relative_error=norm(fe.M-py.M,'fro')/norm(py.M,'fro');
assert(report.stiffness_relative_error<1e-9 && report.mass_relative_error<1e-9,'Element/assembly parity failed.');
report.status='passed';disp(report);
end
