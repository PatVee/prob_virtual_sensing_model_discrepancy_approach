function root=kwmo_setup(root)
% Pass the original package directory; no original files are modified.
assert(isfolder(root),'Original package directory not found.');
[ok,a]=fileattrib(root);assert(ok);root=a.Name;
assert(isfile(fullfile(root,'functions','kwpe_condition.m')),'Wrong original package directory.');
addpath(root,fullfile(root,'functions'));
addpath(genpath(fullfile(root,'dependencies','Stabil')));
addpath(fullfile(root,'dependencies','official_fe'));
if exist('OCTAVE_VERSION','builtin')&&exist('butter','file')~=2,pkg load signal;end
end
