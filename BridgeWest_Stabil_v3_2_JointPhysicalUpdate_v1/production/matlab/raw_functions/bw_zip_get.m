function dest=bw_zip_get(A,name)
i=find(strcmp({A.members.name},name)); assert(numel(i)==1,'Expected one archive member named %s.',name); r=A.members(i);
[~,stem]=fileparts(r.archive); folder=fullfile(A.cache,stem); if ~exist(folder,'dir'), mkdir(folder); end
dest=fullfile(folder,r.name); mark=[dest '.crc']; sig=sprintf('%.0f:%.0f',r.crc,r.bytes);
if exist(dest,'file') && exist(mark,'file'), d=dir(dest); if d.bytes==r.bytes && strcmp(fileread(mark),sig), return; end, end
% Only the selected member is copied; original member paths never become output paths.
z=java.util.zip.ZipFile(r.archive); cleanup=onCleanup(@()z.close()); e=z.getEntry(r.member);
stream=z.getInputStream(e); cleanStream=onCleanup(@()stream.close()); tmp=[dest '.partial'];
f=java.io.File(tmp); if f.exists(), f.delete(); end
options=javaArray('java.nio.file.CopyOption',0);
java.nio.file.Files.copy(stream,f.toPath(),options);
clear cleanStream
info=dir(tmp); assert(info.bytes==r.bytes,'Incomplete extraction: %s',name);
% java.util.zip exposes expected CRC; verify copied bytes before caching.
fid=fopen(tmp,'rb'); assert(fid>=0); c=onCleanup(@()fclose(fid)); crc=java.util.zip.CRC32();
while true
 b=fread(fid,1048576,'*uint8'); if isempty(b), break; end
 crc.update(typecast(b,'int8'));
end
clear c
assert(double(crc.getValue())==r.crc,'ZIP CRC check failed for %s',name);
movefile(tmp,dest,'f'); fid=fopen(mark,'w'); assert(fid>=0); fprintf(fid,'%s',sig); fclose(fid);
end
