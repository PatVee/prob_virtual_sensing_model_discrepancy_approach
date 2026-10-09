function T=bw_read_csv(file)
% Explicit comma delimiter and preserved headers avoid release-dependent
% readtable variable-name changes.
try
    T=readtable(file,'Delimiter',',','VariableNamingRule','preserve','TextType','string');
catch
    T=readtable(file,'Delimiter',',','VariableNamingRule','preserve');
end
end

