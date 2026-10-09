function h=bw_sha256(file)
fid=fopen(file,'rb'); assert(fid>=0,'Cannot read prediction file'); cleanup=onCleanup(@()fclose(fid));
b=fread(fid,Inf,'*uint8'); md=java.security.MessageDigest.getInstance('SHA-256'); md.update(typecast(b,'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
