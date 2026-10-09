function figs=plot_bridgewest_stabil_displacement(outputRoot,varargin)
%PLOT_BRIDGEWEST_STABIL_DISPLACEMENT Plot saved Stabil displacement predictions.
%
% The event figure shows BA and BB vertical-displacement spatial profiles at
% each girder's predicted own-midspan-peak time, with nominal 95% intervals.
% If assessment-only PLS snapshots exist, they are overlaid without changing
% the prediction. A second figure summarises amplitude transfer and spatial
% NRMSE across all available 2025 PLS passages.
%
% Example:
%   plot_bridgewest_stabil_displacement('D:\BW_Stabil_v3_1_results', ...
%       'Event','2025_run40');
assert(nargin>=1&&isfolder(outputRoot),'Result folder does not exist: %s',outputRoot);
ip=inputParser;
addParameter(ip,'Event','2025_run40',@(x)ischar(x)||isstring(x));
addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'Save',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Visible','on',@(x)any(strcmpi(string(x),["on","off"])));
addParameter(ip,'FontSize',13,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(ip,varargin{:});o=ip.Results;o.Event=char(string(o.Event));
if isempty(o.OutputDir),o.OutputDir=fullfile(outputRoot,'figures');end
o.OutputDir=char(string(o.OutputDir));if o.Save&&~isfolder(o.OutputDir),mkdir(o.OutputDir);end
P=load(fullfile(outputRoot,'predictions',[o.Event '.mat']));
figs=struct();
figs.event=figure('Color','w','Visible',o.Visible,'Name',['Displacement ' o.Event], ...
    'Position',[100 100 1280 520]);
tl=tiledlayout(figs.event,1,2,'TileSpacing','compact','Padding','compact');
measFile=fullfile(outputRoot,'displacement_validation',[o.Event '_measured_PLS_snapshot.csv']);
if isfile(measFile),truth=readtable(measFile,'TextType','string');else,truth=table();end
for k=1:2
    g={'BA','BB'};g=g{k};ax=nexttile(tl,k);hold(ax,'on');
    [x,mu,lo,hi,tPeak]=predictedOwnPeak(P,g);co=colororder(ax);c=co(1,:);
    hBand=fill(ax,[x;flipud(x)],[lo;flipud(hi)],c,'FaceAlpha',0.16,'EdgeColor','none', ...
        'DisplayName','Nominal 95% interval');
    hPred=plot(ax,x,mu,'LineWidth',1.8,'DisplayName','Predicted own-peak');
    hMeas=gobjects(0);
    if ~isempty(truth)
        q=truth(string(truth.axis)==string(g),:);
        if ~isempty(q)
            [qx,ord]=sort(double(q.x_m));qm=double(q.measured_mm(ord));
            hMeas=plot(ax,qx,qm,'k.-','LineWidth',1.25,'MarkerSize',8,'DisplayName','Measured PLS own-peak');
        end
    end
    yline(ax,0,'-','Color',[0.65 0.65 0.65],'HandleVisibility','off');
    grid(ax,'on');box(ax,'off');
    xlabel(ax,'Dataset coordinate (m)');ylabel(ax,'Vertical displacement (mm)');
    title(ax,sprintf('(%c) %s: %s, predicted peak t = %.2f s',char(96+k),o.Event,g,tPeak), ...
        'Interpreter','none','FontWeight','normal');
    set(ax,'FontSize',o.FontSize,'TickDir','out');
    if k==1
        if isempty(hMeas),legend(ax,[hBand hPred],{'Nominal 95% interval','Predicted own-peak'},'Box','off','Location','best');
        else,legend(ax,[hBand hPred hMeas],{'Nominal 95% interval','Predicted own-peak','Measured PLS own-peak'},'Box','off','Location','best');end
    end
end
if o.Save,saveFigure(figs.event,fullfile(o.OutputDir,['bw_stabil_displacement_' o.Event]));end

scoreFile=fullfile(outputRoot,'tables','displacement_scores.csv');
if isfile(scoreFile)
    T=readtable(scoreFile,'TextType','string');
    if ~isempty(T)
        figs.summary=figure('Color','w','Visible',o.Visible,'Name','Displacement campaign summary', ...
            'Position',[120 120 1280 520]);
        tl=tiledlayout(figs.summary,1,2,'TileSpacing','compact','Padding','compact');
        metrics={'peak_ratio','spatial_nrmse'};ylab={'Predicted / measured peak magnitude','Spatial NRMSE'};
        ttl={'(a) Displacement amplitude transfer','(b) Displacement profile error'};
        for k=1:2
            ax=nexttile(tl,k);hold(ax,'on');co=colororder(ax);
            for j=1:2
                g={'BA','BB'};g=g{j};q=T(string(T.girder)==string(g),:);if isempty(q),continue;end
                run=zeros(height(q),1);for n=1:height(q),run(n)=eventRun(q.event(n));end
                [run,ord]=sort(run);yy=double(q.(metrics{k})(ord));
                plot(ax,run,yy,'-o','LineWidth',1.1,'MarkerSize',4.5,'Color',co(j,:), ...
                    'MarkerFaceColor',co(j,:),'DisplayName',g);
            end
            if k==1,yline(ax,1,'--','Color',[0.4 0.4 0.4],'HandleVisibility','off');end
            grid(ax,'on');box(ax,'off');xlabel(ax,'2025 passage number');ylabel(ax,ylab{k});
            title(ax,ttl{k},'FontWeight','normal');set(ax,'FontSize',o.FontSize,'TickDir','out');
            if k==1,legend(ax,'Box','off','Location','best');end
        end
        if o.Save,saveFigure(figs.summary,fullfile(o.OutputDir,'bw_stabil_displacement_campaign_summary'));end
    end
end
end

function [x,mu,lo,hi,tPeak]=predictedOwnPeak(P,g)
assert(isfield(P,'displacement_x_dataset_m'), ...
    'Prediction %s has no displacement_x_dataset_m. Re-run predictions with v3.1.',P.event);
x=double(P.displacement_x_dataset_m(:));
Y=double(P.(['disp_' g '_mean_mm']));L=double(P.(['disp_' g '_lo_mm']));U=double(P.(['disp_' g '_hi_mm']));
assert(size(Y,2)==numel(x),'Displacement station count mismatch for %s.',g);
mid=mean([min(x) max(x)]);[~,jm]=min(abs(x-mid));[~,it]=max(abs(Y(:,jm)));t=double(P.time_s(:));tPeak=t(it);
mu=Y(it,:).';lo=L(it,:).';hi=U(it,:).';[x,ord]=sort(x);mu=mu(ord);lo=lo(ord);hi=hi(ord);
end

function r=eventRun(e)
t=regexp(char(string(e)),'run(\d+)$','tokens','once');assert(~isempty(t));r=str2double(t{1});
end

function saveFigure(fig,stem)
set(fig,'ToolBar','none','MenuBar','none');drawnow;
a=findall(fig,'Type','axes');
for k=1:numel(a)
    try,a(k).Toolbar.Visible='off';catch,end
    try,disableDefaultInteractivity(a(k));catch,end
end
exportgraphics(fig,[stem '.png'],'Resolution',300);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
savefig(fig,[stem '.fig']);
end
