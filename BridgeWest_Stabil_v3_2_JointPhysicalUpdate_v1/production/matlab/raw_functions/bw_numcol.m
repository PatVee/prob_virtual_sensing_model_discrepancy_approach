function x=bw_numcol(T,idx)
v=T{:,idx};
if isnumeric(v) || islogical(v)
    x=double(v);
else
    x=str2double(string(v));
end
x=x(:);
end

