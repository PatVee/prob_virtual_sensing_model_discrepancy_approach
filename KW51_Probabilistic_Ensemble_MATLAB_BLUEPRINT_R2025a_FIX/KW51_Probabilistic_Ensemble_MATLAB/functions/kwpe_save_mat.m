function kwpe_save_mat(path,S)
% Write a complete MAT file before atomically publishing its final name.
folder=fileparts(path);if ~exist(folder,'dir'),mkdir(folder);end
tmp=[tempname(folder) '.mat'];save(tmp,'-struct','S','-v7');
[ok,msg]=movefile(tmp,path,'f');assert(ok,'MAT-file save failed: %s',msg);
end
