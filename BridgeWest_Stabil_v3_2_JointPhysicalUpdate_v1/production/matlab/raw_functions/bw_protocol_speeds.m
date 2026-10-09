function s=bw_protocol_speeds(file)
% Robustly read protocol run number and train speed regardless of MATLAB's
% table-variable sanitisation or small metadata header variants.
T=bw_read_csv(file);
iNo=bw_find_col(T,["no","run","runno","runnumber","record","recordno","recordnumber"],"run number");
iSp=bw_find_col(T,["speedkmh","trainspeedkmh","speed","velocitykmh","velocity"],"train speed");
no=round(bw_numcol(T,iNo)); sp=bw_numcol(T,iSp);
good=isfinite(no)&isfinite(sp)&no>=1;
if ~any(good)
    error('BridgeWest:Protocol','No valid run/speed pairs found in %s.',file);
end
nmax=max(100,max(no(good)));
s=nan(1,nmax);
s(no(good))=sp(good);
fprintf('Protocol metadata: %s -> run="%s", speed="%s"\n', ...
    string(file),string(T.Properties.VariableNames{iNo}),string(T.Properties.VariableNames{iSp}));
end

