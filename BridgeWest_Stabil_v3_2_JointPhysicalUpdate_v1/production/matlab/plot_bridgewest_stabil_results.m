function figs=plot_bridgewest_stabil_results(outputRoot,varargin)
%PLOT_BRIDGEWEST_STABIL_RESULTS Strain plots for one or many 2022/2025 runs.
%
% Single passage (the original calling form is retained):
%   plot_bridgewest_stabil_results(outputRoot,'DataRoot',dataRoot,'Event','2025_run22');
%   plot_bridgewest_stabil_results(outputRoot,'Event','2022_run20');
% All saved passages in either campaign or both:
%   plot_bridgewest_stabil_results(outputRoot,'Event','2022','Visible','off');
%   plot_bridgewest_stabil_results(outputRoot,'Event','2025','Visible','off');
%   plot_bridgewest_stabil_results(outputRoot,'Event','all','Visible','off');
% A selected group:
%   plot_bridgewest_stabil_results(outputRoot,'Campaign',2025,'Runs',[20 22 40]);
%   plot_bridgewest_stabil_results(outputRoot,'Event',{'2022_run20','2025_run22'});
%
% Default output is STRAIN. Displacement remains explicitly accessible:
%   ...,'Figures','displacement'
%   ...,'Figures',{'strain','displacement'}
%
% Cached validation_targets files are preferred; DataRoot is optional when
% those exist. No FE calculation, inference, rescaling, or model tuning occurs.
% Non-ok measurements are not scored. ShowReducedQuality=true permits a
% labelled reduced_quality overlay, without turning it into a primary score.
ip=inputParser;
addParameter(ip,'Event',[]);
addParameter(ip,'Campaign','');
addParameter(ip,'Runs',[]);
addParameter(ip,'DataRoot','',@(x)ischar(x)||isstring(x));
addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'Figures','strain');
addParameter(ip,'Save',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Visible','on',@(x)isscalar(string(x))&&any(strcmpi(string(x),["on","off"])));
addParameter(ip,'FontSize',12,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'DPI',300,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=72);
addParameter(ip,'ShowReducedQuality',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'CloseAfterSave',false,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});o=ip.Results;
outputRoot=char(string(outputRoot));o.DataRoot=char(string(o.DataRoot));
here=fileparts(mfilename('fullpath'));addpath(here,fullfile(here,'raw_functions'));
events=bw_select_saved_events(outputRoot,'Event',o.Event,'Campaign',o.Campaign,'Runs',o.Runs);
which=lower(string(o.Figures));which=which(:);
if any(which=="all"),which=["strain";"displacement";"summary"];end
assert(all(ismember(which,["strain","displacement","summary"])), ...
    'BridgeWest:Figures','Figures must be strain, displacement, summary, all, or a list of these.');
if isempty(o.OutputDir),o.OutputDir=fullfile(outputRoot,'figures');end
o.OutputDir=char(string(o.OutputDir));
if o.Save&&~isfolder(o.OutputDir),mkdir(o.OutputDir);end
figs=struct('events',{events},'strain',{{}},'displacement',{{}},'files',{{}});
for k=1:numel(events)
    eid=events{k};
    if any(which=="strain")
        D=bw_stabil_strain_plot_data(outputRoot,eid,o.DataRoot);
        f=makeStrainFigure(D,o);figs.strain{end+1}=f;
        if numel(events)==1,figs.event=f;end
        if o.Save
            stem=fullfile(o.OutputDir,['bw_stabil_strain_' eid]);exportFigure(f,stem,o);
            figs.files{end+1}=[stem '.png'];
        end
        if o.CloseAfterSave&&o.Save,close(f);end
    end
    if any(which=="displacement")
        assert(exist('plot_bridgewest_stabil_displacement','file')==2, ...
            'BridgeWest:MissingPlotter','Install the v3.1/v3.2 displacement plotter.');
        % Assessment-only saved PLS files are used by the displacement plotter.
        % This wrapper does not silently rerun a full-campaign PLS assessment.
        f=plot_bridgewest_stabil_displacement(outputRoot,'Event',eid, ...
            'OutputDir',o.OutputDir,'Save',o.Save,'Visible',o.Visible,'FontSize',o.FontSize);
        figs.displacement{end+1}=f;
        if o.CloseAfterSave&&o.Save
            names=fieldnames(f);for j=1:numel(names),if isgraphics(f.(names{j}),'figure'),close(f.(names{j}));end,end
        end
    end
end
if any(which=="summary")
    sf=fullfile(outputRoot,'tables','target_summary_primary.csv');
    assert(isfile(sf),'BridgeWest:MissingScores','Run score_bridgewest_stabil_results before requesting the summary.');
    S=readtable(sf,'TextType','string');labs=["A058","A070","B031","B081"];Y=nan(4,2);
    for j=1:4
        for k=1:2
            years=[2022 2025];q=double(S.campaign)==years(k)&string(S.target)==labs(j);
            if any(q),assert(nnz(q)==1);Y(j,k)=S.median_nrmse(q);end
        end
    end
    f=figure('Color','w','Visible',o.Visible,'Name','Primary strain scores','ToolBar','none','MenuBar','none');
    ax=axes(f);bar(ax,1:4,Y);set(ax,'XTick',1:4,'XTickLabel',cellstr(labs),'FontSize',o.FontSize);
    ylabel(ax,'Median range-normalised RMSE');xlabel(ax,'Withheld target');grid(ax,'on');
    legend(ax,{'2022','2025'},'Location','best','Box','off');
    title(ax,'Primary strain assessment; non-ok measurements excluded');figs.summary=f;
    if o.Save,exportFigure(f,fullfile(o.OutputDir,'bw_stabil_strain_campaign_summary'),o);end
    if o.CloseAfterSave&&o.Save,close(f);end
end
if o.Save,fprintf('Saved requested plots for %d passage(s) to:\n%s\n',numel(events),o.OutputDir);end
end

function f=makeStrainFigure(D,o)
f=figure('Color','w','Visible',o.Visible,'Name',['Strain ' D.event], ...
    'Position',[80 80 1420 900],'ToolBar','none','MenuBar','none');
tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
header=[D.event ' | A031/B053-conditioned, withheld strain predictions'];
if strcmp(D.event,'2025_run40'),header=[header ' | previously inspected'];end
title(tl,header,'Interpreter','none','FontWeight','normal');
for j=1:4
    ax=nexttile(tl,j);hold(ax,'on');t=D.time_s;
    hp=plot(ax,t,D.mean(:,j),'LineWidth',1.6,'DisplayName','Predicted strain');
    hb=fill(ax,[t;flipud(t)],[D.lo(:,j);flipud(D.hi(:,j))],hp.Color, ...
        'FaceAlpha',0.17,'EdgeColor','none','DisplayName','Nominal 95% interval');
    uistack(hb,'bottom');handles=[hb hp];names={'Nominal 95% interval','Predicted strain'};
    show=D.primary_score(j)||(o.ShowReducedQuality&&D.available(j)&&strcmpi(D.flags{j},'reduced_quality'));
    if show
        hm=plot(ax,t,D.measurement(:,j),'k-','LineWidth',1.15);
        handles=[hb hm hp];names={'Nominal 95% interval','Withheld measurement','Predicted strain'};
        if ~D.primary_score(j),set(hm,'LineStyle',':');names{2}='Reduced-quality measurement (not scored)';end
    end
    if D.primary_score(j)
        label=sprintf('%s | NRMSE %.3f | r %.3f | coverage %.1f%%', ...
            D.labels{j},D.nrmse(j),D.correlation(j),100*D.coverage(j));
    else
        label=sprintf('%s | prediction retained; measurement %s',D.labels{j},D.flags{j});
    end
    title(ax,label,'Interpreter','none','FontWeight','normal');
    if ~D.primary_score(j)
        text(ax,.02,.06,'Not included in primary validation scores','Units','normalized', ...
            'FontSize',max(o.FontSize-2,8),'Interpreter','none');
    end
    xlabel(ax,'Time from bridge entry (s)');ylabel(ax,'Strain (\mu\epsilon)','Interpreter','tex');
    xlim(ax,[t(1) t(end)]);grid(ax,'on');box(ax,'off');set(ax,'FontSize',o.FontSize,'TickDir','out');
    legend(ax,handles,names,'Location','best','Box','off','FontSize',max(o.FontSize-2,8));
    hideToolbar(ax);
end
end
function hideToolbar(ax)
try,ax.Toolbar.Visible='off';catch,end
try,disableDefaultInteractivity(ax);catch,end
end
function exportFigure(f,stem,o)
drawnow;
axs=findall(f,'Type','axes');for j=1:numel(axs),hideToolbar(axs(j));end
exportgraphics(f,[stem '.png'],'Resolution',o.DPI);
exportgraphics(f,[stem '.pdf'],'ContentType','vector');
savefig(f,[stem '.fig']);
end
