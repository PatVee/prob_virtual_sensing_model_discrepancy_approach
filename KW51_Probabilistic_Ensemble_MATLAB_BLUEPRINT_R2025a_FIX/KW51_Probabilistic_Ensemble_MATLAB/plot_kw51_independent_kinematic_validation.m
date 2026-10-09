function [detail,summary]=plot_kw51_independent_kinematic_validation(out,truthFile,varargin)
%PLOT_KW51_INDEPENDENT_KINEMATIC_VALIDATION
% Validate inferred KW51 passage kinematics ONLY against externally supplied,
% event-matched ground truth that was not used anywhere in the inference.
%
% This function deliberately does NOT construct "actual" values from:
%   - sgRAAB / sgRAAL / sgRABB / sgRABL;
%   - rail-energy track selection;
%   - rail-peak timing or the 4.2 m rail-gauge spacing;
%   - Desiro axle geometry;
%   - reference bridge strains;
%   - withheld target strain/displacement measurements;
%   - FE predictions.
%
% Therefore, if no independent truth file is supplied, the function creates
% a template CSV and exits without producing validation plots.
%
% USAGE
%   plot_kw51_independent_kinematic_validation(out)
%   plot_kw51_independent_kinematic_validation(out,truthFile)
%   plot_kw51_independent_kinematic_validation(...,'Visible','on')
%
% REQUIRED/OPTIONAL TRUTH CSV COLUMNS
%   event                       exact traindata_YYYYMMDD_HHMMSS identifier
%   independent_source          true/false (or 1/0); only true rows are used
%   actual_track                A or B; blank when unavailable
%   actual_speed_kmh            numeric; NaN when unavailable
%   actual_acceleration_m_s2    numeric; NaN when unavailable
%   actual_direction            Leuven or Brussels; blank when unavailable
%   source_reference            free text describing the independent source
%
% Direction convention in the FE kinematic vector:
%   +1 = increasing bridge x = Brussels -> Leuven = Leuven-bound
%   -1 = decreasing bridge x = Leuven -> Brussels = Brussels-bound
%
% OUTPUTS
%   detail  : event-level comparison table
%   summary : metric-level validation table
%
% Files are written to:
%   <runFolder>/independent_kinematic_validation/
%
% The function is intentionally conservative: it never substitutes a proxy
% measurement for missing independent truth.

ip=inputParser;
addParameter(ip,'Visible','on',@(x)ischar(x)||isstring(x));
parse(ip,varargin{:});
visible=char(ip.Results.Visible);

runFolder=resolve_run_folder(out);
outFolder=fullfile(runFolder,'independent_kinematic_validation');
if ~exist(outFolder,'dir'),mkdir(outFolder);end

[eventIDs,inferred]=read_inferred_kinematics(runFolder);

if nargin<2 || isempty(truthFile)
    truthFile=fullfile(outFolder,'kw51_independent_kinematics_template.csv');
end
truthFile=char(truthFile);

if exist(truthFile,'file')~=2
    write_truth_template(truthFile,eventIDs);
    fprintf('\nNo independent kinematic truth file was supplied.\n');
    fprintf('Template created:\n  %s\n',truthFile);
    fprintf(['Populate ONLY externally measured/documented event-level truth, ', ...
        'set independent_source=true for those rows, then rerun this function.\n']);
    fprintf('No validation figure has been generated.\n\n');
    detail=table();summary=table();
    return;
end

truth=readtable(truthFile,'TextType','string');
truth=standardise_truth_table(truth);

% Only explicitly declared independent rows are admissible.
truth=truth(truth.independent_source,:);
if isempty(truth)
    fprintf('\nTruth table contains no rows declared independent_source=true.\n');
    fprintf('No validation figure has been generated.\n\n');
    detail=table();summary=table();
    return;
end

% Reject duplicate event truth to avoid ambiguous validation.
[ue,~,ic]=unique(truth.event,'stable');
counts=accumarray(ic,1);
assert(all(counts==1), ...
    'Independent truth contains duplicate event identifiers. Each event must appear once.');

% Exact event matching only.
[tf,loc]=ismember(truth.event,inferred.event);
if any(~tf)
    warning('KW51:IndependentTruthNotInRun', ...
        '%d independent-truth rows do not match this run and will be ignored.',sum(~tf));
end
truth=truth(tf,:);loc=loc(tf);

if isempty(truth)
    fprintf('\nNo independent-truth event matches the selected KW51 run.\n');
    fprintf('No validation figure has been generated.\n\n');
    detail=table();summary=table();
    return;
end

I=inferred(loc,:);

% Build event-level detail table.
detail=table();
detail.event=truth.event;
detail.source_reference=truth.source_reference;
detail.actual_track=truth.actual_track;
detail.inferred_track=I.track;
detail.track_match=valid_string_equal(truth.actual_track,I.track,["A","B"]);
detail.actual_speed_kmh=truth.actual_speed_kmh;
detail.inferred_speed_kmh=I.speed_kmh;
detail.speed_error_kmh=I.speed_kmh-truth.actual_speed_kmh;
detail.actual_acceleration_m_s2=truth.actual_acceleration_m_s2;
detail.inferred_acceleration_m_s2=I.acceleration_m_s2;
detail.acceleration_error_m_s2=I.acceleration_m_s2-truth.actual_acceleration_m_s2;
detail.actual_direction=truth.actual_direction;
detail.inferred_direction=I.direction;
detail.direction_match=valid_string_equal( ...
    truth.actual_direction,I.direction,["Leuven","Brussels"]);

writetable(detail,fullfile(outFolder,'kw_independent_kinematic_validation.csv'));

summaryRows={};

%% Independent speed validation
mask=isfinite(detail.actual_speed_kmh) & isfinite(detail.inferred_speed_kmh);
if any(mask)
    a=detail.actual_speed_kmh(mask);
    p=detail.inferred_speed_kmh(mask);
    stats=numeric_stats(a,p);
    summaryRows(end+1,:)={'speed_kmh',stats.n,stats.rmse,stats.mae,stats.bias,NaN}; %#ok<AGROW>
    make_numeric_validation_figure(a,p, ...
        'Independent measured speed (km/h)', ...
        'Inferred speed (km/h)', ...
        'Independent speed validation',stats, ...
        fullfile(outFolder,'kw_independent_speed_validation'),visible);
else
    fprintf('Independent speed truth: unavailable; speed plot skipped.\n');
end

%% Independent acceleration validation
mask=isfinite(detail.actual_acceleration_m_s2) & ...
    isfinite(detail.inferred_acceleration_m_s2);
if any(mask)
    a=detail.actual_acceleration_m_s2(mask);
    p=detail.inferred_acceleration_m_s2(mask);
    stats=numeric_stats(a,p);
    summaryRows(end+1,:)={'acceleration_m_s2',stats.n,stats.rmse,stats.mae,stats.bias,NaN}; %#ok<AGROW>
    make_numeric_validation_figure(a,p, ...
        'Independent measured acceleration (m/s^2)', ...
        'Inferred acceleration (m/s^2)', ...
        'Independent acceleration validation',stats, ...
        fullfile(outFolder,'kw_independent_acceleration_validation'),visible);
else
    fprintf('Independent acceleration truth: unavailable; acceleration plot skipped.\n');
end

%% Independent track validation
validActual=ismember(upper(strtrim(detail.actual_track)),["A","B"]);
validPred=ismember(upper(strtrim(detail.inferred_track)),["A","B"]);
mask=validActual & validPred;
if any(mask)
    actual=upper(strtrim(detail.actual_track(mask)));
    pred=upper(strtrim(detail.inferred_track(mask)));
    [C,accuracy]=categorical_confusion(actual,pred,["A","B"]);
    summaryRows(end+1,:)={'track',sum(mask),NaN,NaN,NaN,accuracy}; %#ok<AGROW>
    make_confusion_figure(C,["A","B"], ...
        'Independent track validation',accuracy, ...
        fullfile(outFolder,'kw_independent_track_validation'),visible);
else
    fprintf('Independent track truth: unavailable; track plot skipped.\n');
end

%% Independent direction validation
actual=normalise_direction(detail.actual_direction);
pred=normalise_direction(detail.inferred_direction);
mask=ismember(actual,["Leuven","Brussels"]) & ...
    ismember(pred,["Leuven","Brussels"]);
if any(mask)
    [C,accuracy]=categorical_confusion(actual(mask),pred(mask),["Leuven","Brussels"]);
    summaryRows(end+1,:)={'direction',sum(mask),NaN,NaN,NaN,accuracy}; %#ok<AGROW>
    make_confusion_figure(C,["Leuven","Brussels"], ...
        'Independent direction validation',accuracy, ...
        fullfile(outFolder,'kw_independent_direction_validation'),visible);
else
    fprintf('Independent direction truth: unavailable; direction plot skipped.\n');
end

if isempty(summaryRows)
    summary=table();
    fprintf(['\nNo quantity has independently supplied truth. ', ...
        'No validation figures were generated.\n\n']);
    return;
end

summary=cell2table(summaryRows,'VariableNames', ...
    {'metric','n','RMSE','MAE','bias','match_fraction'});
writetable(summary,fullfile(outFolder,'kw_independent_kinematic_validation_summary.csv'));

fprintf('\nIndependent kinematic validation completed.\n');
fprintf('Results folder:\n  %s\n\n',outFolder);
end

%% ========================================================================
function runFolder=resolve_run_folder(out)
if isstruct(out)
    assert(isfield(out,'run_dir'),'The out structure has no run_dir field.');
    runFolder=char(out.run_dir);
elseif ischar(out)||isstring(out)
    runFolder=char(out);
else
    error('First argument must be the run output structure or a run-folder path.');
end
assert(exist(runFolder,'dir')==7,'KW51 run folder does not exist: %s',runFolder);
end

%% ========================================================================
function [eventIDs,T]=read_inferred_kinematics(runFolder)
predFolder=fullfile(runFolder,'predictions');
files=dir(fullfile(predFolder,'traindata_*_prediction.mat'));
assert(~isempty(files),'No prediction MAT files found in %s',predFolder);

names=string({files.name}).';
eventIDs=erase(names,'_prediction.mat');
[eventIDs,order]=sort(eventIDs);
files=files(order);

n=numel(files);
track=strings(n,1);
speed=nan(n,1);
accel=nan(n,1);
direction=strings(n,1);
directionSign=nan(n,1);

for i=1:n
    P=load(fullfile(files(i).folder,files(i).name),'track','kinematics');
    assert(isfield(P,'kinematics')&&numel(P.kinematics)>=4, ...
        'Prediction file %s has no valid kinematics.',files(i).name);
    k=double(P.kinematics(:));
    speed(i)=k(1)*3.6;
    directionSign(i)=sign(k(3));
    accel(i)=k(4);
    if isfield(P,'track'),track(i)=upper(string(P.track));end
    if directionSign(i)>0
        direction(i)="Leuven";
    elseif directionSign(i)<0
        direction(i)="Brussels";
    end
end

T=table(eventIDs,track,speed,accel,directionSign,direction, ...
    'VariableNames',{'event','track','speed_kmh','acceleration_m_s2', ...
    'direction_sign','direction'});
end

%% ========================================================================
function write_truth_template(path,eventIDs)
folder=fileparts(path);if ~isempty(folder)&&exist(folder,'dir')~=7,mkdir(folder);end
n=numel(eventIDs);
event=eventIDs(:);
independent_source=false(n,1);
actual_track=strings(n,1);
actual_speed_kmh=nan(n,1);
actual_acceleration_m_s2=nan(n,1);
actual_direction=strings(n,1);
source_reference=strings(n,1);
T=table(event,independent_source,actual_track,actual_speed_kmh, ...
    actual_acceleration_m_s2,actual_direction,source_reference);
writetable(T,path);
end

%% ========================================================================
function T=standardise_truth_table(T)
vars=string(T.Properties.VariableNames);
required=["event","independent_source","actual_track","actual_speed_kmh", ...
    "actual_acceleration_m_s2","actual_direction","source_reference"];
missing=required(~ismember(required,vars));
assert(isempty(missing), ...
    'Truth CSV is missing required column(s): %s',strjoin(cellstr(missing),', '));

T.event=strtrim(string(T.event));
T.actual_track=upper(strtrim(string(T.actual_track)));
T.actual_direction=normalise_direction(T.actual_direction);
T.source_reference=string(T.source_reference);

% Accept logical, numeric, or common text representations.
x=T.independent_source;
if islogical(x)
    T.independent_source=x;
elseif isnumeric(x)
    T.independent_source=isfinite(x)&x~=0;
else
    s=lower(strtrim(string(x)));
    T.independent_source=ismember(s,["1","true","yes","y"]);
end

T.actual_speed_kmh=to_numeric(T.actual_speed_kmh);
T.actual_acceleration_m_s2=to_numeric(T.actual_acceleration_m_s2);
end

%% ========================================================================
function x=to_numeric(x)
if isnumeric(x),x=double(x);return;end
x=str2double(string(x));
end

%% ========================================================================
function y=normalise_direction(x)
s=lower(strtrim(string(x)));
y=strings(size(s));
y(ismember(s,["leuven","toward leuven","towards leuven","leuven-bound","leuvenbound"]))="Leuven";
y(ismember(s,["brussels","toward brussels","towards brussels","brussels-bound","brusselsbound"]))="Brussels";
end

%% ========================================================================
function m=valid_string_equal(a,b,allowed)
aa=strtrim(string(a));bb=strtrim(string(b));
valid=ismember(aa,allowed)&ismember(bb,allowed);
m=false(size(aa));m(valid)=aa(valid)==bb(valid);
end

%% ========================================================================
function s=numeric_stats(actual,pred)
e=pred-actual;
s=struct();s.n=numel(e);s.rmse=sqrt(mean(e.^2));s.mae=mean(abs(e));s.bias=mean(e);
end

%% ========================================================================
function [C,accuracy]=categorical_confusion(actual,pred,cats)
C=zeros(numel(cats));
for i=1:numel(cats)
    for j=1:numel(cats)
        C(i,j)=sum(actual==cats(i)&pred==cats(j));
    end
end
accuracy=sum(diag(C))/max(sum(C,'all'),1);
end

%% ========================================================================
function make_numeric_validation_figure(actual,pred,xlab,ylab,ttl,stats,base,visible)
fig=figure('Color','w','Visible',visible,'Position',[100 100 620 560]);
ax=axes(fig);hold(ax,'on');box(ax,'on');grid(ax,'on');
scatter(ax,actual,pred,38,'filled');
lo=min([actual(:);pred(:)]);hi=max([actual(:);pred(:)]);
if ~isfinite(lo)||~isfinite(hi),close(fig);return;end
if hi<=lo
    d=max(abs(lo)*0.05,1e-6);lo=lo-d;hi=hi+d;
else
    d=0.05*(hi-lo);lo=lo-d;hi=hi+d;
end
plot(ax,[lo hi],[lo hi],'k--','LineWidth',1.0);
xlim(ax,[lo hi]);ylim(ax,[lo hi]);axis(ax,'square');
xlabel(ax,xlab);ylabel(ax,ylab);title(ax,ttl);
text(ax,0.04,0.96,sprintf('n = %d\nRMSE = %.4g\nMAE = %.4g\nBias = %.4g', ...
    stats.n,stats.rmse,stats.mae,stats.bias), ...
    'Units','normalized','VerticalAlignment','top','BackgroundColor','w','Margin',4);
save_figure(fig,base);close_if_hidden(fig,visible);
end

%% ========================================================================
function make_confusion_figure(C,cats,ttl,accuracy,base,visible)
fig=figure('Color','w','Visible',visible,'Position',[100 100 620 540]);
ax=axes(fig);imagesc(ax,C);axis(ax,'equal');axis(ax,'tight');colorbar(ax);
xticks(ax,1:numel(cats));yticks(ax,1:numel(cats));
xticklabels(ax,cats);yticklabels(ax,cats);
xlabel(ax,'Inferred');ylabel(ax,'Independent actual');
title(ax,sprintf('%s (match = %.1f%%)',ttl,100*accuracy));
mx=max(C(:));if mx==0,mx=1;end
for i=1:size(C,1)
    for j=1:size(C,2)
        text(ax,j,i,sprintf('%d',C(i,j)), ...
            'HorizontalAlignment','center','FontWeight','bold', ...
            'Color',ternary(C(i,j)>0.5*mx,'w','k'));
    end
end
save_figure(fig,base);close_if_hidden(fig,visible);
end

%% ========================================================================
function out=ternary(cond,a,b)
if cond,out=a;else,out=b;end
end

%% ========================================================================
function save_figure(fig,base)
assert(isgraphics(fig,'figure'),'Invalid figure handle.');
axs=findall(fig,'Type','axes');
for k=1:numel(axs)
    if isprop(axs(k),'Toolbar') && ~isempty(axs(k).Toolbar)
        axs(k).Toolbar.Visible='off';
    end
end
drawnow;
exportgraphics(fig,[base '.png'],'Resolution',220,'BackgroundColor','white');
try
    exportgraphics(fig,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
catch
    exportgraphics(fig,[base '.pdf'],'ContentType','image','Resolution',300,'BackgroundColor','white');
end
end

%% ========================================================================
function close_if_hidden(fig,visible)
if strcmpi(visible,'off')&&isgraphics(fig),close(fig);end
end
