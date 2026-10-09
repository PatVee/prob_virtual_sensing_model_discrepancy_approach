function seed=kwpe_seed(event,baseSeed)
bytes=uint8(event);
if exist('OCTAVE_VERSION','builtin')
    h=hash('sha256',char(bytes));
else
    assert(usejava('jvm'),'The standard MATLAB JVM is required for hashing.');
    md=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
    md.update(typecast(bytes(:),'int8'));
    digest=typecast(md.digest(),'uint8');h=lower(reshape(dec2hex(digest,2).',1,[]));
end
seed=mod(double(baseSeed)+hex2dec(h(1:8)),2^32);
end
