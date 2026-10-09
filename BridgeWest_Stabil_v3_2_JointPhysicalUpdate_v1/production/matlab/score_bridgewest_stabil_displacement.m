function [scores,summary]=score_bridgewest_stabil_displacement(dataRoot,outputRoot)
%SCORE_BRIDGEWEST_STABIL_DISPLACEMENT Assessment-only PLS displacement check.
%
% The saved Stabil prediction is never changed. This function opens 2025 PLS
% files only after prediction and compares independent measured and predicted
% own-peak BA/BB spatial profiles. No PLS quantity is used to infer train loads,
% modify shaker-state weights, tune the FE model or change prediction intervals.
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'raw_functions'));
assert(nargin>=2&&isfolder(dataRoot)&&isfolder(outputRoot), ...
    'Call score_bridgewest_stabil_displacement(dataRoot,outputRoot).');
assert(exist('butter','file')==2&&exist('filtfilt','file')==2, ...
    'PLS preprocessing requires Signal Processing Toolbox.');
files=dir(fullfile(outputRoot,'predictions','2025_run*.mat'));
assert(~isempty(files),'No 2025 prediction MAT files found.');
outDir=fullfile(outputRoot,'displacement_validation');
if ~isfolder(outDir),mkdir(outDir);end
cacheDir=fullfile(outDir,'pls_cache');if ~isfolder(cacheDir),mkdir(cacheDir);end
I=indexPLSArchives(dataRoot,cacheDir);
rows=struct('event',{},'girder',{},'n_coordinates',{},'rmse_mm',{}, ...
    'spatial_nrmse',{},'spatial_correlation',{},'peak_ratio',{}, ...
    'predicted_peak_time_s',{},'pls_available',{},'run40_development',{});
status=struct('event',{},'status',{},'reason',{});
for i=1:numel(files)
    P=load(fullfile(files(i).folder,files(i).name));
    eid=char(string(P.event)); run=eventRun(eid);
    try
        plsFile=plsEvent(I,run);
        truth=bw_pls_raw_snapshot(plsFile);
        writetable(truth,fullfile(outDir,[eid '_measured_PLS_snapshot.csv']));
        for g={'BA','BB'}
            [x,mu,lo,hi,tPeak]=predictedOwnPeak(P,g{1});
            q=truth(string(truth.axis)==string(g{1}),:);
            assert(height(q)>=3,'Too few valid %s PLS coordinates for %s.',g{1},eid);
            [qx,ord]=sort(double(q.x_m)); qm=double(q.measured_mm(ord));
            yp=interp1(x,mu,qx,'linear',NaN);
            ylo=interp1(x,lo,qx,'linear',NaN); yhi=interp1(x,hi,qx,'linear',NaN);
            ok=isfinite(qm)&isfinite(yp)&isfinite(ylo)&isfinite(yhi);
            assert(nnz(ok)>=3,'Insufficient overlap between prediction and PLS for %s %s.',eid,g{1});
            qx=qx(ok);qm=qm(ok);yp=yp(ok);ylo=ylo(ok);yhi=yhi(ok);
            rg=max(qm)-min(qm);rmse=sqrt(mean((yp-qm).^2));C=corrcoef(yp,qm);
            if numel(C)>=4,r=C(1,2);else,r=NaN;end
            peakRatio=max(abs(yp))/max(max(abs(qm)),eps);
            rows(end+1)=struct('event',eid,'girder',g{1},'n_coordinates',numel(qm), ... %#ok<AGROW>
                'rmse_mm',rmse,'spatial_nrmse',rmse/max(rg,eps), ...
                'spatial_correlation',r,'peak_ratio',peakRatio, ...
                'predicted_peak_time_s',tPeak,'pls_available',true, ...
                'run40_development',strcmp(eid,'2025_run40'));
            prof=table(qx,qm,yp,ylo,yhi,'VariableNames', ...
                {'x_m','measured_mm','predicted_mm','predicted_lo_mm','predicted_hi_mm'});
            writetable(prof,fullfile(outDir,sprintf('%s_%s_profile.csv',eid,g{1})));
        end
        status(end+1)=struct('event',eid,'status','scored','reason',''); %#ok<AGROW>
    catch ME
        status(end+1)=struct('event',eid,'status','unavailable','reason',ME.message); %#ok<AGROW>
        fprintf(2,'PLS assessment unavailable for %s: %s\n',eid,ME.message);
    end
end
scores=struct2table(rows);
if isempty(scores)
    summary=table();
else
    summary=aggregate(scores);
end
tdir=fullfile(outputRoot,'tables');if ~isfolder(tdir),mkdir(tdir);end
writetable(scores,fullfile(tdir,'displacement_scores.csv'));
writetable(summary,fullfile(tdir,'displacement_summary.csv'));
writeJson(fullfile(outDir,'displacement_status.json'),status);
meta=struct('prediction_stage_uses_PLS',false,'PLS_first_opened_after_prediction',true, ...
    'comparison','independent measured/predicted own-midspan-peak spatial profiles', ...
    'n_scored_profiles',height(scores),'n_2025_prediction_files',numel(files));
writeJson(fullfile(outDir,'displacement_assessment_summary.json'),meta);
fprintf('Displacement assessment complete: %d BA/BB profiles scored.\n',height(scores));
if ~isempty(summary),disp(summary);end
end

function [x,mu,lo,hi,tPeak]=predictedOwnPeak(P,g)
assert(isfield(P,'displacement_x_dataset_m'),'Prediction %s has no displacement coordinate vector. Re-run with v3.1.',P.event);
x=double(P.displacement_x_dataset_m(:));
meanField=['disp_' g '_mean_mm'];loField=['disp_' g '_lo_mm'];hiField=['disp_' g '_hi_mm'];
assert(isfield(P,meanField)&&isfield(P,loField)&&isfield(P,hiField),'Prediction %s lacks %s displacement fields.',P.event,g);
Y=double(P.(meanField));L=double(P.(loField));U=double(P.(hiField));
assert(size(Y,2)==numel(x),'Displacement station count mismatch in %s %s.',P.event,g);
mid=mean([min(x) max(x)]);[~,jm]=min(abs(x-mid));
[~,it]=max(abs(Y(:,jm)));t=double(P.time_s(:));tPeak=t(it);
mu=Y(it,:).';lo=L(it,:).';hi=U(it,:).';
[x,ord]=sort(x);mu=mu(ord);lo=lo(ord);hi=hi(ord);
end

function I=indexPLSArchives(dataRoot,cacheRoot)
D=dir(fullfile(dataRoot,'**','ltw24_*.zip'));
assert(~isempty(D),'No ltw24_*.zip archives found below %s.',dataRoot);
I=struct('cache',cacheRoot,'members',struct('archive',{},'member',{},'name',{},'bytes',{},'crc',{}));
for i=1:numel(D)
    zpath=fullfile(D(i).folder,D(i).name);z=java.util.zip.ZipFile(zpath);cleanup=onCleanup(@()z.close());en=z.entries();
    while en.hasMoreElements()
        e=en.nextElement();if e.isDirectory(),continue;end
        member=char(e.getName());[~,base,ext]=fileparts(strrep(member,'\','/'));name=[base ext];
        if ~strcmpi(ext,'.h5')||~startsWith(string(name),'PLSW_ICETD_'),continue;end
        row=struct('archive',zpath,'member',member,'name',name,'bytes',double(e.getSize()),'crc',double(e.getCrc()));
        if isempty(I.members)||~any(strcmp({I.members.name},name)),I.members(end+1)=row;end %#ok<AGROW>
    end
    clear cleanup
end
assert(~isempty(I.members),'No PLSW HDF5 members found in LT-W24 archives.');
end

function file=plsEvent(I,run)
pat=sprintf('^PLSW_ICETD_%02d_DT_2025_.*[.]h5$',run);
hit=find(~cellfun(@isempty,regexpi({I.members.name},pat)));
assert(numel(hit)==1,'Expected one PLS record for 2025 run %02d; found %d.',run,numel(hit));
r=I.members(hit);[~,stem]=fileparts(r.archive);folder=fullfile(I.cache,stem);if ~isfolder(folder),mkdir(folder);end
file=fullfile(folder,r.name);mark=[file '.crc'];sig=sprintf('%.0f:%.0f',r.crc,r.bytes);
if isfile(file)&&isfile(mark)
    d=dir(file);if d.bytes==r.bytes&&strcmp(fileread(mark),sig),return;end
end
z=java.util.zip.ZipFile(r.archive);cleanup=onCleanup(@()z.close());e=z.getEntry(r.member);stream=z.getInputStream(e);cleanStream=onCleanup(@()stream.close());
tmp=[file '.partial'];if isfile(tmp),delete(tmp);end
options=javaArray('java.nio.file.CopyOption',0);java.nio.file.Files.copy(stream,java.io.File(tmp).toPath(),options);clear cleanStream
info=dir(tmp);assert(info.bytes==r.bytes,'Incomplete PLS extraction for %s.',r.name);movefile(tmp,file,'f');
fid=fopen(mark,'w');assert(fid>=0);fprintf(fid,'%s',sig);fclose(fid);clear cleanup
end

function r=eventRun(eid)
t=regexp(eid,'run(\d+)$','tokens','once');assert(~isempty(t),'Cannot parse run from %s.',eid);r=str2double(t{1});
end

function S=aggregate(T)
out=struct('girder',{},'n_passages',{},'median_spatial_nrmse',{}, ...
    'median_spatial_correlation',{},'median_peak_ratio',{},'median_rmse_mm',{});
for g={'BA','BB'}
    q=strcmp(T.girder,g{1});if ~any(q),continue;end
    out(end+1)=struct('girder',g{1},'n_passages',sum(q), ... %#ok<AGROW>
        'median_spatial_nrmse',median(T.spatial_nrmse(q),'omitnan'), ...
        'median_spatial_correlation',median(T.spatial_correlation(q),'omitnan'), ...
        'median_peak_ratio',median(T.peak_ratio(q),'omitnan'), ...
        'median_rmse_mm',median(T.rmse_mm(q),'omitnan'));
end
S=struct2table(out);
end

function writeJson(file,x)
fid=fopen(file,'w');assert(fid>=0,'Cannot write %s.',file);c=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(x,'PrettyPrint',true));
end
