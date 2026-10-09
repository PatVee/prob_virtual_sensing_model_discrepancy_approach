function kw51_write_csv(path,headers,rows)
fid=fopen(path,'w');assert(fid>=0,'Cannot write %s',path);clean=onCleanup(@()fclose(fid));
allrows=[headers(:).';rows];
for i=1:size(allrows,1)
 for j=1:size(allrows,2)
  v=allrows{i,j};if isnumeric(v)||islogical(v),s=sprintf('%.16g',double(v));else,s=char(v);end
  s=strrep(s,'"','""');if j>1,fprintf(fid,',');end;fprintf(fid,'"%s"',s);
 end
 fprintf(fid,'\n');
end
end
