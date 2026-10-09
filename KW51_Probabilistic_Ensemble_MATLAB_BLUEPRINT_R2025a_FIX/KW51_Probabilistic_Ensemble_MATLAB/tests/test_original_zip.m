function test_original_zip(zipFile)
% Optional end-to-end smoke test on the ORIGINAL December part-1 ZIP archive.
% Generates its own inference/validation/figures; does not use saved results.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'functions'));
test_kw51_probabilistic();
output=tempname;mkdir(output);
out=run_kw51_probabilistic_matlab(zipFile,output, ...
    'Events',{'traindata_20191210_170156'},'CandidateCount',90, ...
    'PairedCount',12,'PosteriorDrawCount',64,'Validate',true,'Plot',true);
V=load(fullfile(out.run_dir,'validation','traindata_20191210_170156_validation.mat'));
assert(all(isfinite(V.filtered_mean(:))),'Invalid prediction.');
assert(all(V.upper95(:)>=V.lower95(:)),'Reversed uncertainty bounds.');
F=load(fullfile(root,'tests','fe_fixture.mat'));
E=load(fullfile(out.run_dir,'inputs','traindata_20191210_170156_input.mat'));
assert(max(abs(E.y_ref(:)-F.E_full.y_ref(:)))<1e-6,'Original-ZIP preprocessing mismatch.');
assert(max(abs(E.rail_offsets_s(:)-F.E_full.rail_offsets_s(:)))<1e-12,'Rail-offset preprocessing mismatch.');
fprintf('Original-ZIP smoke test completed: %s\n',out.run_dir);
fprintf('This smoke test uses reduced ensemble sizes, not the full 900/48/256 settings.\n');
end
