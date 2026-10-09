function bw_write_json(path,value)
fid=fopen(path,'w'); assert(fid>=0,'Cannot write %s',path); cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value));
end
