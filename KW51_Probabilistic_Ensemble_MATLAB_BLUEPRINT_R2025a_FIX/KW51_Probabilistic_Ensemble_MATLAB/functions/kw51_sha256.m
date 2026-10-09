function h=kw51_sha256(path)
fid=fopen(path,'rb');assert(fid>=0,'Cannot read %s',path);clean=onCleanup(@()fclose(fid));
if exist('OCTAVE_VERSION','builtin')
 bytes=fread(fid,Inf,'*uint8');h=hash('sha256',char(bytes.'));return
end
assert(usejava('jvm'),'KW51 requires MATLAB with its standard JVM for SHA-256 verification.');
md=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
while ~feof(fid)
 bytes=fread(fid,1024*1024,'*uint8');if ~isempty(bytes),md.update(typecast(bytes(:),'int8'));end
end
bytes=typecast(md.digest(),'uint8');h=lower(reshape(dec2hex(bytes,2).',1,[]));
end
