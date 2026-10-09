function [scores,summary,refs]=score_bridgewest_stabil_results(dataRoot,outputRoot)
%SCORE_BRIDGEWEST_STABIL_RESULTS Separate protected-target assessment stage.
%
% This function is deliberately separate from the inference runner. It opens
% A058/A070/B031/B081 only after prediction files already exist. Their values
% cannot alter loads, FE state, shaker discrepancy or prediction intervals.
%
% All four target predictions are retained for every predicted passage. A
% reduced-quality measurement can still receive an assessment row, but only
% quality_flag == 'ok' contributes to the primary aggregate summary.
root=fileparts(mfilename('fullpath'));pkg=fileparts(fileparts(root));
addpath(root,fullfile(root,'raw_functions'));
assert(nargin>=2&&~isempty(dataRoot)&&~isempty(outputRoot),'Call score_bridgewest_stabil_results(dataRoot,outputRoot).');
c=jsondecode(fileread(fullfile(pkg,'production','config','prediction_config.json')));
pre=jsondecode(fileread(fullfile(pkg,'production','config','preprocessing_config.json')));
A=bw_prediction_zip_index(dataRoot,fullfile(outputRoot,'validation_archive_cache'));
files=dir(fullfile(outputRoot,'predictions','*_run*.mat'));assert(~isempty(files),'No prediction MAT files found.');

rows=struct('event',{},'campaign',{},'target',{},'n',{},'rmse_microstrain',{},'nrmse',{}, ...
    'correlation',{},'coverage',{},'range_ratio',{},'mean_bias_microstrain',{}, ...
    'quality_flag',{},'assessment_available',{},'primary_score',{},'run40_development',{});
rrows=struct('event',{},'campaign',{},'reference',{},'nrmse',{},'correlation',{},'coverage',{});
modelTargets=double(c.target_model_indices(:)).';targetLabels=cellstr(string(c.target_labels));
measuredTargets=cellstr(string(c.measured_target_channels));
modelRefs=double(c.reference_model_indices(:)).';measuredRefs=cellstr(string(c.measured_reference_channels));

for i=1:numel(files)
    P=load(fullfile(files(i).folder,files(i).name));
    campaign=str2double(P.event(1:4));
    assert(isfield(P,'entry_s'),'Prediction %s predates strict v3 and has no frozen entry time.',P.event);

    % Protected targets are first read here, after the prediction exists.
    Z=bw_prepare_validation_targets(A,P.event,pre,outputRoot,measuredTargets, ...
        P.speed_m_s,P.entry_s,P.time_s);
    chs=cellstr(string(Z.channels));flags=cellstr(string(Z.flags));

    for j=1:numel(modelTargets)
        col=find(strcmp(chs,measuredTargets{j}),1);
        assert(~isempty(col),'Validation target %s was not prepared for %s.',measuredTargets{j},P.event);
        flag=flags{col};yy=Z.Y(:,col);
        mm=P.strain_mean_microstrain(:,modelTargets(j));
        lo=P.strain_lo_microstrain(:,modelTargets(j));hi=P.strain_hi_microstrain(:,modelTargets(j));
        available=all(isfinite(yy))&&numel(yy)==numel(mm);
        if available,rr=responseScore(yy,mm,lo,hi);else,rr=emptyScore();end
        primary=available&&strcmp(flag,'ok');
        rows(end+1)=struct('event',P.event,'campaign',campaign,'target',targetLabels{j}, ... %#ok<AGROW>
            'n',numel(yy),'rmse_microstrain',rr.rmse,'nrmse',rr.nrmse,'correlation',rr.r, ...
            'coverage',rr.coverage,'range_ratio',rr.range_ratio,'mean_bias_microstrain',rr.bias, ...
            'quality_flag',flag,'assessment_available',available,'primary_score',primary, ...
            'run40_development',strcmp(P.event,'2025_run40'));
    end

    % Reference goodness-of-fit is computed from the permitted-input file,
    % not by reopening target-containing traffic arrays.
    rf=fullfile(outputRoot,'permitted_inputs',[P.event '.mat']);assert(exist(rf,'file')==2,'Missing permitted input %s.',P.event);
    R=load(rf);auditReferenceOnly(R,measuredRefs,measuredTargets,P.event);
    for j=1:numel(modelRefs)
        yy=R.reference_strain_microstrain(:,j);
        mm=P.strain_mean_microstrain(:,modelRefs(j));
        lo=P.strain_lo_microstrain(:,modelRefs(j));hi=P.strain_hi_microstrain(:,modelRefs(j));
        rr=responseScore(yy,mm,lo,hi);
        rrows(end+1)=struct('event',P.event,'campaign',campaign,'reference',measuredRefs{j}, ... %#ok<AGROW>
            'nrmse',rr.nrmse,'correlation',rr.r,'coverage',rr.coverage);
    end
end

scores=struct2table(rows);refs=struct2table(rrows);
assert(~isempty(scores),'No target assessment rows were produced.');
expectedRows=4*numel(files);assert(height(scores)==expectedRows,'Expected four target rows per predicted passage.');
primary=scores(scores.primary_score,:);
summary=aggregate(primary);
allAvailable=scores(scores.assessment_available,:);summaryAll=aggregate(allAvailable);

% A separate view excludes run40 because it was previously inspected. The
% principal requested 36-passage tables still include it and label it.
uninspected=primary(~primary.run40_development,:);summaryUninspected=aggregate(uninspected);

tdir=fullfile(outputRoot,'tables');if ~exist(tdir,'dir'),mkdir(tdir);end
writetable(scores,fullfile(tdir,'target_assessment_all.csv'));
writetable(primary,fullfile(tdir,'target_scores_primary.csv'));
writetable(summary,fullfile(tdir,'target_summary_primary.csv'));
writetable(summaryAll,fullfile(tdir,'target_summary_all_available.csv'));
writetable(summaryUninspected,fullfile(tdir,'target_summary_uninspected_only.csv'));
writetable(refs,fullfile(tdir,'reference_scores.csv'));

assessmentSummary=struct('predicted_passages',numel(files),'target_predictions',expectedRows, ...
    'target_measurements_available',height(allAvailable),'primary_quality_qualified_rows',height(primary), ...
    'reference_channels',{measuredRefs},'withheld_targets',{measuredTargets}, ...
    'targets_read_only_after_prediction',true, ...
    'run40_included_in_requested_36_passage_tables',true, ...
    'run40_also_flagged_as_previously_inspected',true);
bw_write_json(fullfile(outputRoot,'assessment_summary.json'),assessmentSummary);

fprintf('Assessment complete: %d predictions across %d passages.\n',expectedRows,numel(files));
fprintf('%d rows are quality-qualified for the primary aggregate.\n',height(primary));
disp(summary);
end

function auditReferenceOnly(R,refs,targets,eid)
assert(isfield(R,'reference_names')&&isfield(R,'reference_strain_microstrain'),'Permitted input %s is incomplete.',eid);
names=cellstr(string(R.reference_names));
assert(isequal(names(:),refs(:)),'Permitted input %s does not contain exactly A031/B053 in the configured order.',eid);
assert(size(R.reference_strain_microstrain,2)==numel(refs),'Permitted input %s has unexpected structural channel count.',eid);
for j=1:numel(targets)
    assert(~any(strcmp(names,targets{j})),'Protected target %s leaked into permitted inputs for %s.',targets{j},eid);
end
if isfield(R,'contains_protected_target_measurements')
    assert(~logical(R.contains_protected_target_measurements),'Permitted input %s declares protected target content.',eid);
end
end

function r=responseScore(y,m,lo,hi)
y=double(y(:));m=double(m(:));lo=double(lo(:));hi=double(hi(:));
rg=max(y)-min(y);rmse=sqrt(mean((m-y).^2));C=corrcoef(y,m);
if numel(C)>=4,cr=C(1,2);else,cr=NaN;end
r=struct('rmse',rmse,'nrmse',rmse/max(rg,eps),'r',cr,'coverage',mean(y>=lo&y<=hi), ...
    'range_ratio',(max(m)-min(m))/max(rg,eps),'bias',mean(m-y));
end

function r=emptyScore()
r=struct('rmse',NaN,'nrmse',NaN,'r',NaN,'coverage',NaN,'range_ratio',NaN,'bias',NaN);
end

function S=aggregate(T)
out=struct('campaign',{},'target',{},'n_passages',{},'median_nrmse',{}, ...
    'median_correlation',{},'mean_coverage',{},'median_range_ratio',{});
if isempty(T),S=struct2table(out);return;end
camps=unique(T.campaign);targets=unique(string(T.target),'stable');
for ic=1:numel(camps)
    for it=1:numel(targets)
        q=T.campaign==camps(ic)&string(T.target)==targets(it);if ~any(q),continue;end
        out(end+1)=struct('campaign',camps(ic),'target',char(targets(it)), ... %#ok<AGROW>
            'n_passages',sum(q),'median_nrmse',median(T.nrmse(q),'omitnan'), ...
            'median_correlation',median(T.correlation(q),'omitnan'), ...
            'mean_coverage',mean(T.coverage(q),'omitnan'), ...
            'median_range_ratio',median(T.range_ratio(q),'omitnan'));
    end
end
S=struct2table(out);
end
