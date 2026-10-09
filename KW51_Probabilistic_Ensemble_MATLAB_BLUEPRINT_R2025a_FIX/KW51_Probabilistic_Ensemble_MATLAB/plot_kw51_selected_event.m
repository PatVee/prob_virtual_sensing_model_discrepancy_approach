function files = plot_kw51_selected_event(out,eventSelector,varargin)
%PLOT_KW51_SELECTED_EVENT Plot strain and displacement for one KW51 passage.
%
% files = plot_kw51_selected_event(out,eventSelector)
%
% INPUTS
%   out           : output struct returned by run_kw51_probabilistic_matlab, or
%                   the run-folder path.
%   eventSelector : preferred form is a chronological accepted-event number
%                   1,2,...,N. A full event ID is also accepted, e.g.
%                   'traindata_20191201_171118'.
%
% Event numbers are assigned after sorting the accepted event IDs
% chronologically. Therefore event 1 is the earliest accepted passage in
% the current run, event 2 is the next accepted passage, and so on.
%
% NAME/VALUE
%   'Visible'      : 'on' or 'off' (default 'on').
%   'OutputFolder' : default <runFolder>/selected_event_figures.
%   'OutputName'   : default 'kw_selected_<eventID>'.
%
% The function reads saved validation products only. It does not refit the
% FE model, loads, discrepancy, or posterior weights.

ip=inputParser;
addParameter(ip,'Visible','on');
addParameter(ip,'OutputFolder','');
addParameter(ip,'OutputName','');
parse(ip,varargin{:});
opt=ip.Results;

runFolder=resolve_run_folder(out);
[eventID,eventNumber,eventCount]=resolve_event_selector(runFolder,eventSelector);
if isempty(opt.OutputFolder)
    outputFolder=fullfile(runFolder,'selected_event_figures');
else
    outputFolder=char(opt.OutputFolder);
end
if ~exist(outputFolder,'dir'),mkdir(outputFolder);end

settingsFile=fullfile(runFolder,'settings.mat');
validationFile=fullfile(runFolder,'validation',[eventID '_validation.mat']);
resultsFile=fullfile(runFolder,'validation_results.mat');
assert(exist(settingsFile,'file')==2,'Missing settings.mat in %s.',runFolder);
assert(exist(validationFile,'file')==2, ...
    'No saved validation file for event %s.',eventID);

c=load(settingsFile);
V=load(validationFile);
if exist(resultsFile,'file')==2
    R=load(resultsFile,'scores');
    scores=R.scores;
else
    scores={};
end

targets=as_cellstr(c.protected_names);
t=V.time_s(:);
idxDefault=find(abs(V.translation_sigma_values-c.translation_sd_mm)<1e-12,1);
assert(~isempty(idxDefault),'Default translation scale is absent from validation file.');

blue=[0.08 0.37 0.51];
pale=[0.75 0.85 0.91];

fig=figure('Color','w','Visible',opt.Visible,'Units','centimeters', ...
    'Position',[2 2 29.7 18.2]);
layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

for j=1:4
    ax=nexttile(layout,j);hold(ax,'on');

    fill(ax,[t;flipud(t)], ...
        [V.lower95(:,j,idxDefault);flipud(V.upper95(:,j,idxDefault))], ...
        pale,'EdgeColor','none','FaceAlpha',0.65, ...
        'DisplayName','Nominal 95% interval');

    plot(ax,t,V.filtered_mean(:,j),'-','Color',blue,'LineWidth',1.45, ...
        'DisplayName','Posterior mean');

    if iscell(V.native_time)
        tx=V.native_time{j}; yy=V.native_truth{j};
    else
        tx=[];yy=[];
    end
    if ~isempty(tx)
        if j<=2
            plot(ax,tx,yy,'k--','LineWidth',1.0, ...
                'DisplayName','Withheld measurement');
        else
            plot(ax,tx,yy,'k.--','LineWidth',0.9,'MarkerSize',8, ...
                'DisplayName','Withheld measurement');
        end
    end

    startTime=V.kinematics(2)-1;
    finishTime=V.kinematics(2)+116.316/V.kinematics(1)+max(V.rail_offsets_s)+2;
    xlim(ax,[startTime finishTime]);
    grid(ax,'on');box(ax,'off');
    set(ax,'FontName','Times New Roman','FontSize',10,'LineWidth',0.75);
    title(ax,sprintf('(%c) %s',char('a'+j-1),targets{j}), ...
        'Interpreter','none','FontWeight','normal');
    xlabel(ax,'Time from record start (s)');
    if j<=2
        ylabel(ax,'Strain (\muepsilon)');
    else
        ylabel(ax,'Displacement (mm)');
    end

    [nr,coverage]=event_score(scores,V.event_id,targets{j},c.translation_sd_mm);
    if isfinite(nr)
        text(ax,0.02,0.97,sprintf('NRMSE = %.3f; coverage = %.1f%%',nr,100*coverage), ...
            'Units','normalized','VerticalAlignment','top','FontName','Times New Roman', ...
            'FontSize',8.5,'BackgroundColor','w','Margin',2);
    end

    if j==3
        legend(ax,'Location','best','FontName','Times New Roman','FontSize',8);
    end
end

if isfinite(eventNumber)
    heading=sprintf('Event %d/%d: %s | Track %s | KW51 blueprint BC', ...
        eventNumber,eventCount,strrep(V.event_id,'_','\_'),V.track);
else
    heading=sprintf('%s | Track %s | KW51 blueprint BC', ...
        strrep(V.event_id,'_','\_'),V.track);
end
title(layout,heading, ...
    'FontName','Times New Roman','FontSize',11,'FontWeight','normal');

if isempty(opt.OutputName)
    if isfinite(eventNumber)
        outputName=sprintf('kw_selected_event_%02d_%s',eventNumber,eventID);
    else
        outputName=['kw_selected_' eventID];
    end
else
    outputName=char(opt.OutputName);
end
base=fullfile(outputFolder,outputName);
save_figure(fig,base);

files=struct('pdf',[base '.pdf'],'png',[base '.png'],'figure',fig);
if strcmpi(opt.Visible,'off'),close(fig);end
fprintf('Selected-event figure written to:\n  %s\n  %s\n',files.pdf,files.png);
end

function runFolder=resolve_run_folder(out)
if isstruct(out)
    assert(isfield(out,'run_dir'),'Input struct has no run_dir field.');
    runFolder=char(out.run_dir);
else
    runFolder=char(out);
end
assert(exist(runFolder,'dir')==7,'Run folder does not exist: %s',runFolder);
end

function [eventID,eventNumber,eventCount]=resolve_event_selector(runFolder,eventSelector)
catalogFile=fullfile(runFolder,'input_catalog.mat');
assert(exist(catalogFile,'file')==2,'Missing input_catalog.mat in %s.',runFolder);
C=load(catalogFile,'inputs');
assert(isfield(C,'inputs')&&~isempty(C.inputs),'No accepted events found in input_catalog.mat.');
names={C.inputs.event_id};
names=sort(names); % IDs contain YYYYMMDD_HHMMSS, so lexical order is chronological.
eventCount=numel(names);

if isnumeric(eventSelector)
    assert(isscalar(eventSelector)&&isfinite(eventSelector)&&eventSelector==round(eventSelector), ...
        'Numeric event selector must be one integer in the range 1:%d.',eventCount);
    eventNumber=double(eventSelector);
    assert(eventNumber>=1&&eventNumber<=eventCount, ...
        'Event number must be between 1 and %d.',eventCount);
    eventID=names{eventNumber};
elseif ischar(eventSelector)||(isstring(eventSelector)&&isscalar(eventSelector))
    eventID=char(eventSelector);
    eventNumber=find(strcmp(names,eventID),1);
    assert(~isempty(eventNumber),'Event ID %s is not an accepted event in this run.',eventID);
else
    error('eventSelector must be an integer event number or a full event ID.');
end

fprintf('Selected event %d of %d: %s\n',eventNumber,eventCount,eventID);
end

function [nr,coverage]=event_score(scores,eventID,target,sigma)
nr=NaN;coverage=NaN;
if isempty(scores),return;end
try
    ix=strcmp(scores(:,1),eventID) & strcmp(scores(:,4),target) & ...
        abs(cell2mat(scores(:,6))-sigma)<1e-12;
    k=find(ix,1);
    if ~isempty(k)
        nr=scores{k,8};coverage=scores{k,11};
    end
catch
    nr=NaN;coverage=NaN;
end
end

function c=as_cellstr(x)
if iscell(x),c=x(:);elseif isstring(x),c=cellstr(x(:));else,c=cellstr(x);end
end

function save_figure(fig,base)
% Robust MATLAB R2025a-compatible figure export.
assert(isgraphics(fig,'figure'),'Invalid or deleted figure handle passed to save_figure.');

% Force pending graphics operations to complete before export.
drawnow;

% Hide interactive axes toolbars so they are not embedded in exports.
axs=findall(fig,'Type','axes');
for k=1:numel(axs)
    try
        if isprop(axs(k),'Toolbar') && ~isempty(axs(k).Toolbar)
            axs(k).Toolbar.Visible='off';
        end
    catch
        % Toolbar availability varies across MATLAB graphics versions.
    end
end

drawnow;

pngFile=[base '.png'];
pdfFile=[base '.pdf'];

% exportgraphics is preferred over print for figures containing tiledlayout
% objects in current MATLAB releases (including R2025a).
try
    exportgraphics(fig,pngFile, ...
        'Resolution',220, ...
        'BackgroundColor','white');
catch MEpng
    warning('KW51:PNGExportFallback', ...
        'exportgraphics PNG failed (%s). Falling back to getframe/imwrite.',MEpng.message);
    frame=getframe(fig);
    imwrite(frame.cdata,pngFile);
end

try
    exportgraphics(fig,pdfFile, ...
        'ContentType','vector', ...
        'BackgroundColor','white');
catch MEvec
    warning('KW51:PDFVectorFallback', ...
        'Vector PDF export failed (%s). Retrying as image PDF.',MEvec.message);
    exportgraphics(fig,pdfFile, ...
        'ContentType','image', ...
        'Resolution',220, ...
        'BackgroundColor','white');
end
end