function summary=validate_kw51_probabilistic_matlab(runFolder)
% Reads original protected observations only after verifying frozen inference.
% Complete conditional paths are filtered/interpolated before forming bounds.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'functions'));
if isstruct(runFolder),runFolder=runFolder.run_dir;end
kwpe_verify_freeze(runFolder);
I=load(fullfile(runFolder,'input_catalog.mat'));inputs=I.inputs;
c=load(fullfile(runFolder,'settings.mat'));targets=c.protected_names(:).';
assert(abs(c.translation_sd_mm)<1e-12&&numel(c.translation_sensitivity_mm)==1&& ...
    abs(c.translation_sensitivity_mm(1))<1e-12, ...
    'This validation routine expects the blueprint BC with no common-translation prior.');
sigmas=0;L=1;
folder=fullfile(runFolder,'validation');if ~exist(folder,'dir'),mkdir(folder);end
scores={};skipped={};
for k=1:numel(inputs)
    P=load(fullfile(runFolder,'predictions',[inputs(k).event_id '_prediction.mat']));
    t=P.time_s(:);fsPred=1/median(diff(t));baselinePred=t>=0&t<1;
    assert(any(baselinePred),'No prediction baseline samples.');
    raw=load(inputs(k).raw_file,'predat_sg','predat_d');
    V=struct('event_id',P.event_id,'split',P.split,'track',P.track,'time_s',t, ...
        'labels',{targets},'kinematics',P.kinematics,'rail_offsets_s',P.rail_offsets_s, ...
        'translation_sigma_values',sigmas,'display_translation_sd_mm',c.translation_sd_mm, ...
        'posterior_ess',P.posterior_ess);
    V.truth=nan(numel(t),4);V.filtered_mean=V.truth;
    V.lower95=nan(numel(t),4,L);V.upper95=V.lower95;
    V.native_time=cell(1,4);V.native_truth=cell(1,4);V.native_prediction=cell(1,4);
    V.native_lower95=cell(4,L);V.native_upper95=cell(4,L);V.native_score_mask=cell(1,4);
    V.validation_cutoff_Hz=nan(1,4);
    for j=1:4
        if j<=2,container='predat_sg';else,container='predat_d';end
        if ~isfield(raw,container),skipped(end+1,:)={P.event_id,targets{j},'Missing measurement container'};continue;end %#ok<AGROW>
        a=raw.(container);names=cellstr(a.labels);index=find(strcmp(names,targets{j}));
        if numel(index)~=1,skipped(end+1,:)={P.event_id,targets{j},'Missing or duplicate channel'};continue;end %#ok<AGROW>
        y=double(a.tdata(:,index));fs=double(a.fs);tt=(double(a.sdn(:))-P.strain_start_sdn)*86400;
        if any(~isfinite(y)),skipped(end+1,:)={P.event_id,targets{j},'Nonfinite validation data'};continue;end %#ok<AGROW>
        if j<=2,cutoff=10;scale=1e6;units='microstrain';else,cutoff=min(.8,.4*fs);scale=1e3;units='mm';end
        baseline=tt>=0&tt<1;
        if ~any(baseline),skipped(end+1,:)={P.event_id,targets{j},'No predeclared baseline samples'};continue;end %#ok<AGROW>
        y=kw51_filter(y,cutoff,fs)*scale;y=y-median(y(baseline));
        V.validation_cutoff_Hz(j)=cutoff;V.truth(:,j)=interp1(tt,y,t,'linear',NaN);
        pg=kw51_filter(P.mean(:,j+3),cutoff,fsPred);pg=pg-median(pg(baselinePred));
        V.filtered_mean(:,j)=pg;
        % Stored paths are S x T x 4; filtering expects T x S.
        paths=reshape(double(P.raw_target_draws(:,:,j)),size(P.raw_target_draws,1),numel(t)).';
        paths=kw51_filter(paths,cutoff,fsPred);paths=paths-median(paths(baselinePred,:),1);
        % Blueprint BC: Brussels provides the longitudinal reference.
        % No independent common-translation process is added to bearings.
        if j<=2,ts=t;yy=V.truth(:,j);
        else,valid=tt>=t(1)&tt<=t(end);ts=tt(valid);yy=y(valid);end
        pp=interp1(t,pg,ts);
        V.native_time{j}=ts;V.native_truth{j}=yy;V.native_prediction{j}=pp;
        finish=P.kinematics(2)+116.316/P.kinematics(1)+max(P.rail_offsets_s)+2;
        good=ts>=P.kinematics(2)-1&ts<=finish&isfinite(yy)&isfinite(pp);
        V.native_score_mask{j}=good;
        if sum(good)<8||max(yy(good))-min(yy(good))<1e-9
            skipped(end+1,:)={P.event_id,targets{j},'Insufficient valid range/samples'};continue %#ok<AGROW>
        end
        truth=yy(good);prediction=pp(good);range=max(truth)-min(truth);
        rmse=sqrt(mean((prediction-truth).^2));nr=rmse/range;
        corr=corrcoef(prediction,truth);r=corr(1,2);
        nse=1-sum((prediction-truth).^2)/sum((truth-mean(truth)).^2);
        peak=100*abs(max(abs(prediction))-max(abs(truth)))/max(abs(truth));
        for h=1:L
            draw=paths;
            sd=std(draw,0,2);V.lower95(:,j,h)=pg-c.interval_z*sd;V.upper95(:,j,h)=pg+c.interval_z*sd;
            % Interpolate each path, not the marginal SD or interval endpoints.
            nativeDraw=interp1(t,draw,ts);nativeSD=std(nativeDraw,0,2);
            lower=pp-c.interval_z*nativeSD;upper=pp+c.interval_z*nativeSD;
            V.native_lower95{j,h}=lower;V.native_upper95{j,h}=upper;
            lo=lower(good);hi=upper(good);
            coverage=mean(truth>=lo&truth<=hi);width=mean(hi-lo);
            intervalScore=mean(hi-lo+40*max(lo-truth,0)+40*max(truth-hi,0));
            scores(end+1,:)={P.event_id,P.split,P.track,targets{j},sum(good),sigmas(h), ...
                rmse,nr,r,nse,coverage,width,width/range,intervalScore/range,peak,P.posterior_ess,units}; %#ok<AGROW>
        end
    end
    kwpe_save_mat(fullfile(folder,[P.event_id '_validation.mat']),V);
    fprintf('Validated %d/%d: %s\n',k,numel(inputs),P.event_id);
end
headers={'event','split','track','target','n_samples','translation_sd_mm','RMSE','NRMSE', ...
    'r','NSE','coverage95','mean_interval_width','normalized_interval_width','normalized_interval_score', ...
    'peak_amplitude_error_percent','ESS','units'};
kw51_write_csv(fullfile(runFolder,'validation_scores.csv'),headers,scores);
kw51_write_csv(fullfile(runFolder,'validation_skipped.csv'),{'event','target','reason'},skipped);
summary={};splits={'development','evaluation'};
for s=1:numel(splits)
    for j=1:4
        for h=1:L
            if isempty(scores),continue;end
            ix=strcmp(scores(:,2),splits{s})&strcmp(scores(:,4),targets{j})&cell2mat(scores(:,6))==sigmas(h);
            if ~any(ix),continue;end
            nr=cell2mat(scores(ix,8));rr=cell2mat(scores(ix,9));cover=cell2mat(scores(ix,11));
            width=cell2mat(scores(ix,13));properScore=cell2mat(scores(ix,14));ess=cell2mat(scores(ix,16));
            summary(end+1,:)={splits{s},targets{j},sigmas(h),numel(nr),median(nr),mean(nr), ...
                median(rr),mean(cover),median(width),median(properScore),median(ess),sum(nr<=.15)}; %#ok<AGROW>
            if sigmas(h)==c.translation_sd_mm
                fprintf('%s %s: median NRMSE %.4f; mean coverage %.1f%%; width/range %.2f\n', ...
                    splits{s},targets{j},median(nr),100*mean(cover),median(width));
            end
        end
    end
end
kw51_write_csv(fullfile(runFolder,'validation_summary.csv'), ...
    {'split','target','translation_sd_mm','n','median_NRMSE','mean_NRMSE','median_r', ...
    'mean_coverage95','median_normalized_width','median_interval_score','median_ESS','n_NRMSE_le_0_15'},summary);
kwpe_save_mat(fullfile(runFolder,'validation_results.mat'), ...
    struct('score_headers',{headers},'scores',{scores},'summary',{summary},'skipped',{skipped}));
kwpe_verify_freeze(runFolder);
end
