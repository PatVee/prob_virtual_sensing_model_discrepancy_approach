function figureFolder=plot_kw51_probabilistic_matlab(out,varargin)
% Separate plotting only. Reads saved validation/prediction files; no refit.
% Default: first evaluation passage on each track, plus validation summary.
% To plot every passage: plot_kw51_probabilistic_matlab(out,'Events','all')
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'functions'));
if isstruct(out),runFolder=out.run_dir;else,runFolder=char(out);end
ip=inputParser;addParameter(ip,'Events','examples');addParameter(ip,'Visible','off');parse(ip,varargin{:});opt=ip.Results;
assert(exist(fullfile(runFolder,'validation_results.mat'),'file')==2,'Run validation before plotting measured comparisons.');
c=load(fullfile(runFolder,'settings.mat'));R=load(fullfile(runFolder,'validation_results.mat'));
scores=R.scores;summary=R.summary;targets=c.protected_names(:).';
figureFolder=fullfile(runFolder,'figures');if ~exist(figureFolder,'dir'),mkdir(figureFolder);end
files=kw51_files(fullfile(runFolder,'validation'),'_validation\.mat$');files=sort(files);
selection=opt.Events;if isstring(selection),selection=cellstr(selection);end
if ischar(selection),selection={selection};end
chosen={};seenTracks={};
for k=1:numel(files)
    V=load(files{k},'event_id','split','track');
    if any(strcmp(selection,'all')),chosen{end+1}=files{k}; %#ok<AGROW>
    elseif any(strcmp(selection,'examples'))
        if strcmp(V.split,'evaluation')&&~any(strcmp(seenTracks,V.track))
            chosen{end+1}=files{k};seenTracks{end+1}=V.track; %#ok<AGROW>
        end
    elseif any(strcmp(selection,V.event_id)),chosen{end+1}=files{k}; %#ok<AGROW>
    end
end
if isempty(chosen)&&~isempty(files)
    warning('KW51:PlotSubset','No evaluation example found; showing the first available passage.');chosen=files(1);
end
blue=[.08 .37 .51];pale=[.75 .85 .91];
for k=1:numel(chosen)
    V=load(chosen{k});t=V.time_s(:);
    h=find(abs(V.translation_sigma_values-c.translation_sd_mm)<1e-12,1);
    fig=figure('Color','w','Visible',opt.Visible,'Position',[80 80 1200 760]);
    for j=1:4
        ax=subplot(2,2,j);hold(ax,'on');
        fill(ax,[t;flipud(t)],[V.lower95(:,j,h);flipud(V.upper95(:,j,h))],pale, ...
            'EdgeColor','none','DisplayName','Nominal 95% interval');
        plot(ax,t,V.filtered_mean(:,j),'Color',blue,'LineWidth',1.5,'DisplayName','Ensemble prediction');
        if ~isempty(V.native_time{j})
            if j>=3,lineStyle='k.--';else,lineStyle='k--';end
            plot(ax,V.native_time{j},V.native_truth{j},lineStyle,'LineWidth',1.0,'MarkerSize',8,'DisplayName','Withheld measurement');
        end
        start=V.kinematics(2)-1;finish=V.kinematics(2)+116.316/V.kinematics(1)+max(V.rail_offsets_s)+2;
        xlim(ax,[start finish]);grid(ax,'on');box(ax,'off');title(ax,targets{j},'Interpreter','none');
        xlabel(ax,'Time from record start (s)');
        if j<=2,ylabel(ax,'Strain (microstrain)');else,ylabel(ax,'Displacement (mm)');end
        if ~isempty(scores)
            ix=strcmp(scores(:,1),V.event_id)&strcmp(scores(:,4),targets{j})& ...
                abs(cell2mat(scores(:,6))-c.translation_sd_mm)<1e-12;
            a=find(ix,1);if ~isempty(a)
                text(ax,.02,.97,sprintf('NRMSE %.3f | coverage %.1f%%',scores{a,8},100*scores{a,11}), ...
                    'Units','normalized','VerticalAlignment','top','BackgroundColor','w','FontSize',9);
            end
        end
        if j==3,legend(ax,'Location','best','FontSize',8);end
    end
    sgtitle(sprintf('KW51 %s | Track %s | blueprint BC',V.event_id,V.track),'Interpreter','none');
    kwpe_save_figure(fig,fullfile(figureFolder,[V.event_id '_ensemble']));close(fig);
end
% Summaries show all evaluation passages and do not select by error.
if isempty(summary),return;end
split='evaluation';if ~any(strcmp(summary(:,1),split)),split='development';end
values=nan(4,3);
for j=1:4
    ix=strcmp(summary(:,1),split)&strcmp(summary(:,2),targets{j})& ...
        abs(cell2mat(summary(:,3))-c.translation_sd_mm)<1e-12;
    a=find(ix,1);if ~isempty(a),values(j,:)=[summary{a,5},100*summary{a,8},summary{a,9}];end
end
fig=figure('Color','w','Visible',opt.Visible,'Position',[80 80 1200 430]);
titles={'Median prediction error','Nominal interval coverage','Interval width'};
ylabels={'Range-normalised RMSE','Mean coverage (%)','Median width / measured range'};
for j=1:3
    ax=subplot(1,3,j);bar(ax,values(:,j),'FaceColor',blue);hold(ax,'on');
    if j==1,plot(ax,[.5 4.5],[.15 .15],'k:');elseif j==2,plot(ax,[.5 4.5],[95 95],'k--');ylim(ax,[0 110]);end
    set(ax,'XTick',1:4,'XTickLabel',targets,'TickLabelInterpreter','none','XTickLabelRotation',25);
    ylabel(ax,ylabels{j});title(ax,titles{j});grid(ax,'on');box(ax,'off');
end
sgtitle(sprintf('KW51 %s set - coverage must be interpreted together with interval width',split));
kwpe_save_figure(fig,fullfile(figureFolder,'KW51_ensemble_validation_summary'));close(fig);
fprintf('PNG and vector PDF figures: %s\n',figureFolder);
end
