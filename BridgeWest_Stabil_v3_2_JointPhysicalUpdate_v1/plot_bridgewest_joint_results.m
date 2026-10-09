function figs = plot_bridgewest_joint_results(packageRoot,v3ResultRoot,outputRoot,varargin)
%PLOT_BRIDGEWEST_JOINT_RESULTS Plot completed joint-discrepancy Bridge West results.
%
% This function DOES NOT rerun FE analysis or inference. It reads only the
% frozen outputs created by run_bridgewest_joint_calibrated_v3_2 and the
% assessment-only validation files in BridgeWest_Stabil_Shaker_v3_results.
%
% Example:
%   figs = plot_bridgewest_joint_results( ...
%       packageRoot,resultRoot,outputRoot,'Event','2025_run40');
%
% Optional name/value arguments:
%   'Event'      selected event for detailed plots (default 2025_run40);
%                use 'all' to plot every completed 2025 passage
%   'OutputDir'  figure directory (default outputRoot/figures_joint)
%   'Visible'    'on' or 'off' (default 'on')
%   'Save'       true/false (default true)
%   'DPI'        PNG resolution (default 300)
%   'FontSize'   base font size (default 12)
%
% Figures produced:
%   1. 2025 total inferred load and independent static total
%   2. selected-event bogie load comparison
%   3. selected-event withheld strain histories
%   4. selected-event BA/BB displacement profiles
%   5. campaign target-strain NRMSE
%   6. campaign displacement peak ratio and spatial NRMSE
%   7. 2022 frozen strain-recovery calibration

ip=inputParser;
addParameter(ip,'Event','2025_run40',@(x)ischar(x)||isstring(x));
addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'Visible','on',@(x)any(strcmpi(string(x),["on","off"])));
addParameter(ip,'Save',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'DPI',300,@(x)isnumeric(x)&&isscalar(x)&&x>=72);
addParameter(ip,'FontSize',12,@(x)isnumeric(x)&&isscalar(x)&&x>0);
% Internal option used when 'Event','all' recursively generates event figures.
addParameter(ip,'SummaryPlots',true,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});o=ip.Results;

packageRoot=char(string(packageRoot));
v3ResultRoot=char(string(v3ResultRoot));
outputRoot=char(string(outputRoot));
eid=char(string(o.Event));
assert(isfolder(packageRoot),'packageRoot not found: %s',packageRoot);
assert(isfolder(v3ResultRoot),'resultRoot not found: %s',v3ResultRoot);
assert(isfolder(outputRoot),'outputRoot not found: %s',outputRoot);
if isempty(o.OutputDir),o.OutputDir=fullfile(outputRoot,'figures_joint');end
o.OutputDir=char(string(o.OutputDir));
if o.Save && ~isfolder(o.OutputDir),mkdir(o.OutputDir);end

cfgFile=fullfile(packageRoot,'production','config','prediction_config.json');
assert(isfile(cfgFile),'Missing prediction config: %s',cfgFile);
c=jsondecode(fileread(cfgFile));

feFile=fullfile(outputRoot,'fe_model','bridgewest_lorenzen_fe.mat');
assert(isfile(feFile),'Missing completed FE file: %s',feFile);
F=load(feFile,'fe');fe=F.fe;
xDisp=double(fe.observation_maps.x_dataset_m(:));

loadFile=fullfile(outputRoot,'2025_load_validation.csv');
strainFile=fullfile(outputRoot,'2025_strain_validation.csv');
dispFile=fullfile(outputRoot,'2025_displacement_validation.csv');
calFile=fullfile(outputRoot,'2022_reference_scale_by_event.csv');
assert(isfile(loadFile),'Missing %s. Complete the joint run first.',loadFile);
L=readtable(loadFile,'TextType','string');
if isfile(strainFile),ST=readtable(strainFile,'TextType','string');else,ST=table();end
if isfile(dispFile),DT=readtable(dispFile,'TextType','string');else,DT=table();end
if isfile(calFile),CT=readtable(calFile,'TextType','string');else,CT=table();end

figs=struct();

% 'all' means: make campaign summaries once, then event-level plots for every
% completed 2025 prediction. No FE analysis or inference is rerun.
if strcmpi(strtrim(eid),'all')
    pf=dir(fullfile(outputRoot,'predictions','2025_run*_prediction.mat'));
    assert(~isempty(pf),'No completed 2025 prediction files found in %s',fullfile(outputRoot,'predictions'));
    events=erase(string({pf.name}),'_prediction.mat');
    runsAll=arrayfun(@(x)eventRun(x),events);
    [~,io]=sort(runsAll);events=events(io);
    figs=struct('events',struct());
    for ie=1:numel(events)
        makeSummary = (ie==1) && o.SummaryPlots;
        fThis=plot_bridgewest_joint_results(packageRoot,v3ResultRoot,outputRoot, ...
            'Event',events(ie),'OutputDir',o.OutputDir,'Visible',o.Visible, ...
            'Save',o.Save,'DPI',o.DPI,'FontSize',o.FontSize,'SummaryPlots',makeSummary);
        figs.events.(matlab.lang.makeValidName(char(events(ie))))=fThis;
    end
    if o.Save
        fprintf('Saved detailed figures for %d completed 2025 passages to:\n%s\n',numel(events),o.OutputDir);
    end
    return;
end

%% 1. Campaign total load
if o.SummaryPlots
figs.load_campaign=figure('Color','w','Visible',o.Visible,'Name','Joint discrepancy: 2025 load totals', ...
    'Position',[100 100 1120 520],'ToolBar','none','MenuBar','none');
ax=axes(figs.load_campaign);hold(ax,'on');
runs=cellfun(@eventRun,cellstr(string(L.event)));
[runs,ord]=sort(runs);Q=double(L.total_kN(ord));qlo=double(L.total_q025_kN(ord));qhi=double(L.total_q975_kN(ord));
static=double(L.static_total_kN(ord));
for i=1:numel(runs)
    plot(ax,[runs(i) runs(i)],[qlo(i) qhi(i)],'-','LineWidth',1.15,'HandleVisibility','off');
end
h1=plot(ax,runs,Q,'o-','LineWidth',1.35,'MarkerSize',5,'DisplayName','Reference-inferred total');
h2=plot(ax,runs,static,'k--','LineWidth',1.35,'DisplayName','Independent static total');
grid(ax,'on');box(ax,'off');xlabel(ax,'2025 passage number');ylabel(ax,'Total bogie load (kN)');
title(ax,'Joint discrepancy: inferred train load with 95% posterior total interval','FontWeight','normal');
legend(ax,[h1 h2],'Location','best','Box','off');styleAxes(ax,o);
if o.Save,saveFigure(figs.load_campaign,fullfile(o.OutputDir,'joint_2025_total_load'),o);end
end

%% Load selected event prediction
predFile=fullfile(outputRoot,'predictions',[eid '_prediction.mat']);
assert(isfile(predFile),'Selected event prediction not found: %s',predFile);
S=load(predFile,'pred','provenance');p=S.pred;
pinFile=fullfile(v3ResultRoot,'permitted_inputs',[eid '.mat']);assert(isfile(pinFile),'Missing %s',pinFile);
pin=load(pinFile);t=double(pin.time_s(:));

%% 2. Selected-event bogie load
figs.event_load=figure('Color','w','Visible',o.Visible,'Name',['Joint load ' eid], ...
    'Position',[120 120 980 520],'ToolBar','none','MenuBar','none');
ax=axes(figs.event_load);hold(ax,'on');
q=double(p.load_mean_kN(:));sd=double(p.load_sd_kN(:));z=double(c.interval_z);
lo=max(double(c.load_lower_kN),q-z*sd);hi=min(double(c.load_upper_kN),q+z*sd);
W=readtable(fullfile(packageRoot,'audit','verified_static_bogie_loads.csv'));
year=str2double(eid(1:4));Ws=sortrows(W(W.campaign==year,:),'bogie');staticB=double(Ws.static_load_kN(:));
for b=1:8,plot(ax,[b b],[lo(b) hi(b)],'-','LineWidth',1.3,'HandleVisibility','off');end
h1=plot(ax,1:8,q,'o-','LineWidth',1.35,'MarkerSize',6,'DisplayName','Posterior mean');
h2=plot(ax,1:8,staticB,'ks--','LineWidth',1.15,'MarkerSize',5,'DisplayName','Campaign static weighing');
xlim(ax,[0.5 8.5]);xticks(ax,1:8);xlabel(ax,'Bogie number');ylabel(ax,'Bogie resultant (kN)');
grid(ax,'on');box(ax,'off');title(ax,[eid ' — bogie-load inference'],'Interpreter','none','FontWeight','normal');
legend(ax,[h1 h2],'Location','best','Box','off');styleAxes(ax,o);
text(ax,0.02,0.04,'Vertical ranges are marginal mean \pm 1.96 SD, clipped to prior bounds.', ...
    'Units','normalized','FontSize',max(o.FontSize-2,8));
if o.Save,saveFigure(figs.event_load,fullfile(o.OutputDir,['joint_bogie_load_' eid]),o);end

%% 3. Selected-event withheld strain histories
figs.event_strain=figure('Color','w','Visible',o.Visible,'Name',['Joint strain ' eid], ...
    'Position',[80 80 1350 820],'ToolBar','none','MenuBar','none');
tl=tiledlayout(figs.event_strain,2,2,'TileSpacing','compact','Padding','compact');
title(tl,[eid ' — withheld strain; A031/B053-conditioned joint discrepancy'], ...
    'Interpreter','none','FontWeight','normal');
Vfile=fullfile(v3ResultRoot,'validation_targets',[eid '.mat']);
if isfile(Vfile),V=load(Vfile);else,V=[];end
for j=1:4
    ax=nexttile(tl,j);hold(ax,'on');midx=double(c.target_model_indices(j));
    mu=double(p.strain_mean(:,midx));plo=double(p.strain_lo(:,midx));phi=double(p.strain_hi(:,midx));
    fill(ax,[t;flipud(t)],[plo;flipud(phi)],[0.75 0.82 0.90],'EdgeColor','none','FaceAlpha',0.45, ...
        'DisplayName','Nominal 95% interval');
    if ~isempty(V) && size(V.Y,2)>=j
        yy=double(V.Y(:,j));plot(ax,t,yy,'k','LineWidth',1.05,'DisplayName','Measured withheld strain');
    end
    plot(ax,t,mu,'LineWidth',1.4,'DisplayName','Joint prediction');
    grid(ax,'on');box(ax,'off');xlabel(ax,'Time from bridge entry (s)');ylabel(ax,'Strain (\mu\epsilon)');
    title(ax,sprintf('(%c) %s',char(96+j),char(string(c.target_labels(j)))),'FontWeight','normal');styleAxes(ax,o);
    if j==1,legend(ax,'Location','best','Box','off');end
end
if o.Save,saveFigure(figs.event_strain,fullfile(o.OutputDir,['joint_strain_' eid]),o);end

%% 4. Selected-event BA/BB displacement profiles
figs.event_displacement=figure('Color','w','Visible',o.Visible,'Name',['Joint displacement ' eid], ...
    'Position',[100 100 1280 540],'ToolBar','none','MenuBar','none');
tl=tiledlayout(figs.event_displacement,1,2,'TileSpacing','compact','Padding','compact');
plsFile=fullfile(v3ResultRoot,'displacement_validation',[eid '_measured_PLS_snapshot.csv']);
if isfile(plsFile),PLS=readtable(plsFile,'TextType','string');else,PLS=table();end
for ig=1:2
    ax=nexttile(tl,ig);hold(ax,'on');g={'BA','BB'};g=g{ig};
    Y=double(p.(['disp_' g '_mean']));Lo=double(p.(['disp_' g '_lo']));Hi=double(p.(['disp_' g '_hi']));
    mid=mean([min(xDisp) max(xDisp)]);[~,jm]=min(abs(xDisp-mid));[~,it]=max(abs(Y(:,jm)));
    [xs,io]=sort(xDisp);mu=Y(it,io).';plo=Lo(it,io).';phi=Hi(it,io).';
    fill(ax,[xs;flipud(xs)],[plo;flipud(phi)],[0.75 0.82 0.90],'EdgeColor','none','FaceAlpha',0.45, ...
        'DisplayName','Nominal 95% interval');
    plot(ax,xs,mu,'LineWidth',1.5,'DisplayName','Predicted own-peak');
    if ~isempty(PLS)
        qq=PLS(strcmpi(string(PLS.axis),g),:);
        if ~isempty(qq)
            if ismember('x_m',qq.Properties.VariableNames),xm=double(qq.x_m);ym=double(qq.measured_mm);else,xm=double(qq.x);ym=double(qq.measured);end
            [xm,oo]=sort(xm);ym=ym(oo);plot(ax,xm,ym,'k.-','LineWidth',1.15,'MarkerSize',7,'DisplayName','Measured PLS own-peak');
        end
    end
    yline(ax,0,'-','Color',[0.65 0.65 0.65],'HandleVisibility','off');grid(ax,'on');box(ax,'off');
    xlabel(ax,'Dataset coordinate (m)');ylabel(ax,'Vertical displacement (mm)');
    title(ax,sprintf('(%c) %s; predicted peak t = %.2f s',char(96+ig),g,t(it)),'FontWeight','normal');styleAxes(ax,o);
    if ig==1,legend(ax,'Location','best','Box','off');end
end
if o.Save,saveFigure(figs.event_displacement,fullfile(o.OutputDir,['joint_displacement_' eid]),o);end

%% 5--7. Campaign summaries and frozen 2022 calibration
if o.SummaryPlots
%% 5. Campaign strain NRMSE
if ~isempty(ST)
    figs.strain_summary=figure('Color','w','Visible',o.Visible,'Name','Joint strain summary', ...
        'Position',[120 120 1120 520],'ToolBar','none','MenuBar','none');
    ax=axes(figs.strain_summary);hold(ax,'on');targets=string(c.target_labels(:));
    for j=1:numel(targets)
        qq=ST(string(ST.target)==targets(j),:);if isempty(qq),continue;end
        rr=cellfun(@eventRun,cellstr(string(qq.event)));[rr,io]=sort(rr);plot(ax,rr,double(qq.nrmse(io)),'-o','LineWidth',1.1,'MarkerSize',4.5,'DisplayName',targets(j));
    end
    grid(ax,'on');box(ax,'off');xlabel(ax,'2025 passage number');ylabel(ax,'Range-normalised RMSE');
    title(ax,'Withheld-strain error across the 2025 campaign','FontWeight','normal');legend(ax,'Location','best','Box','off');styleAxes(ax,o);
    if o.Save,saveFigure(figs.strain_summary,fullfile(o.OutputDir,'joint_2025_strain_nrmse'),o);end
end

%% 6. Campaign displacement summary
if ~isempty(DT)
    figs.displacement_summary=figure('Color','w','Visible',o.Visible,'Name','Joint displacement summary', ...
        'Position',[120 120 1280 520],'ToolBar','none','MenuBar','none');
    tl=tiledlayout(figs.displacement_summary,1,2,'TileSpacing','compact','Padding','compact');
    metrics={'peak_ratio','spatial_nrmse'};ylab={'Predicted / measured peak magnitude','Spatial NRMSE'};
    ttl={'(a) Displacement amplitude transfer','(b) Displacement profile error'};
    for k=1:2
        ax=nexttile(tl,k);hold(ax,'on');
        for ig=1:2
            g={'BA','BB'};g=g{ig};qq=DT(strcmpi(string(DT.girder),g),:);if isempty(qq),continue;end
            rr=cellfun(@eventRun,cellstr(string(qq.event)));[rr,io]=sort(rr);plot(ax,rr,double(qq.(metrics{k})(io)),'-o','LineWidth',1.1,'MarkerSize',4.5,'DisplayName',g);
        end
        if k==1,yline(ax,1,'k--','HandleVisibility','off');end
        grid(ax,'on');box(ax,'off');xlabel(ax,'2025 passage number');ylabel(ax,ylab{k});title(ax,ttl{k},'FontWeight','normal');styleAxes(ax,o);
        if k==1,legend(ax,'Location','best','Box','off');end
    end
    if o.Save,saveFigure(figs.displacement_summary,fullfile(o.OutputDir,'joint_2025_displacement_summary'),o);end
end

%% 7. Frozen 2022 strain-recovery scale
if ~isempty(CT)
    figs.calibration=figure('Color','w','Visible',o.Visible,'Name','2022 strain-recovery calibration', ...
        'Position',[140 140 1060 500],'ToolBar','none','MenuBar','none');
    ax=axes(figs.calibration);hold(ax,'on');rr=cellfun(@eventRun,cellstr(string(CT.event)));[rr,io]=sort(rr);ss=double(CT.scale_event(io));
    plot(ax,rr,ss,'o-','LineWidth',1.2,'MarkerSize',5,'DisplayName','2022 event scale');
    med=median(ss,'omitnan');yline(ax,med,'k--','LineWidth',1.25,'DisplayName',sprintf('Frozen median = %.4f',med));
    grid(ax,'on');box(ax,'off');xlabel(ax,'2022 passage number');ylabel(ax,'Strain-recovery scale');
    title(ax,'A031/B053 + independent 2022 static weighing calibration','FontWeight','normal');legend(ax,'Location','best','Box','off');styleAxes(ax,o);
    if o.Save,saveFigure(figs.calibration,fullfile(o.OutputDir,'joint_2022_strain_recovery_calibration'),o);end
end
end % o.SummaryPlots

if o.Save,fprintf('Saved joint-discrepancy figures to:\n%s\n',o.OutputDir);end
end

function r=eventRun(e)
tok=regexp(char(string(e)),'run(\d+)$','tokens','once');assert(~isempty(tok),'Cannot parse run from %s',string(e));r=str2double(tok{1});
end
function styleAxes(ax,o)
set(ax,'FontSize',o.FontSize,'TickDir','out','LineWidth',0.8);try,ax.Toolbar.Visible='off';catch,end
end
function saveFigure(fig,stem,o)
drawnow;exportgraphics(fig,[stem '.png'],'Resolution',o.DPI);exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig(fig,[stem '.fig']);
end
