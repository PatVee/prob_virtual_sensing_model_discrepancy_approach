function files=kw51_files(folder,pattern)
% Portable recursive listing (also works without dir('**')).
files={};items=dir(folder);
for k=1:numel(items)
 n=items(k).name;if strcmp(n,'.')||strcmp(n,'..')||strcmp(n,'__MACOSX'),continue;end
 p=fullfile(folder,n);
 if items(k).isdir,files=[files;kw51_files(p,pattern)]; %#ok<AGROW>
 elseif ~isempty(regexp(n,pattern,'once')),files{end+1,1}=p;end %#ok<AGROW>
end
end
