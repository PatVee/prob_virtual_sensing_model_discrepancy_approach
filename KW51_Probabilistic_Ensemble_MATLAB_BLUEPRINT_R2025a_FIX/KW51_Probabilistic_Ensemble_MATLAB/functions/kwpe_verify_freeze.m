function kwpe_verify_freeze(runFolder)
F=load(fullfile(runFolder,'prediction_freeze.mat'));
for k=1:numel(F.files)
    assert(exist(F.files{k},'file')==2,'Missing frozen file: %s',F.files{k});
    assert(strcmp(F.hashes{k},kw51_sha256(F.files{k})),'KW51:FrozenChanged', ...
        'Frozen file changed: %s',F.files{k});
end
end
