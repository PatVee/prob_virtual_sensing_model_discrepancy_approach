function figureFolder = plot_kw51_manuscript_figures(out,varargin)
%PLOT_KW51_MANUSCRIPT_FIGURES Rebuild all KW51 manuscript figures in MATLAB.
%
% figureFolder = plot_kw51_manuscript_figures(out)
%
% INPUT
%   out : output struct returned by run_kw51_probabilistic_matlab, or the
%         saved run-folder path.
%
% NAME/VALUE
%   'Visible'      : 'on' or 'off' (default 'off').
%   'OutputFolder' : default <runFolder>/manuscript_figures.
%
% FIGURES CREATED
%   kw_fe_mesh.pdf/png
%   kw_inputs.pdf/png
%   kw_posterior_load_benchmark.pdf/png
%   kw_posterior_load_benchmark_17passage_summary.pdf/png
%   kw_ensemble_validation_summary.pdf/png
%   kw_prediction_histories_trackB.pdf/png
%   kw_prediction_histories_trackA.pdf/png
%
% This is plotting/post-processing only. It reads the frozen outputs from
% the probabilistic run and does not refit the bridge or use withheld data
% to alter inference.

ip=inputParser;
addParameter(ip,'Visible','off');
addParameter(ip,'OutputFolder','');
parse(ip,varargin{:});
opt=ip.Results;

runFolder=resolve_run_folder(out);
sourceRoot=package_root();
addpath(sourceRoot,fullfile(sourceRoot,'functions'));

% Use the exact run-local FE bank that generated the saved draws. Never
% replay posterior weights with a different boundary condition.
bankFile=fullfile(runFolder,'fe_bank','model_15.mat');
assert(exist(bankFile,'file')==2, ...
    'Missing run-local FE bank: %s. Re-run inference with the blueprint package.',bankFile);
Bcheck=load(bankFile,'G','boundary_description');
assert(size(Bcheck.G,2)==2, ...
    'Run-local FE bank is not the two-coordinate blueprint bank.');
fprintf('Manuscript plotting FE bank: %s\n',bankFile);
fprintf('Support discrepancy coordinates: %d\n',size(Bcheck.G,2));
if isfield(Bcheck,'boundary_description'),fprintf('BC: %s\n',Bcheck.boundary_description);end

if isempty(opt.OutputFolder)
    figureFolder=fullfile(runFolder,'manuscript_figures');
else
    figureFolder=char(opt.OutputFolder);
end
if ~exist(figureFolder,'dir'),mkdir(figureFolder);end

assert(exist(fullfile(runFolder,'validation_results.mat'),'file')==2, ...
    'Validation results are required. Run with Validate=true first.');

c=load(fullfile(runFolder,'settings.mat'));
R=load(fullfile(runFolder,'validation_results.mat'));
I=load(fullfile(runFolder,'input_catalog.mat'));
inputs=I.inputs;

make_fe_mesh(sourceRoot,figureFolder,opt.Visible);
make_inputs_figure(runFolder,inputs,figureFolder,opt.Visible);
loadSummary=make_load_benchmark_figures(runFolder,bankFile,inputs,figureFolder,opt.Visible);
make_validation_summary(R.summary,c,figureFolder,opt.Visible);

% Exact two passage histories used in the manuscript.
plot_kw51_selected_event(runFolder,'traindata_20191201_080952', ...
    'Visible',opt.Visible,'OutputFolder',figureFolder, ...
    'OutputName','kw_prediction_histories_trackB');
plot_kw51_selected_event(runFolder,'traindata_20191201_171118', ...
    'Visible',opt.Visible,'OutputFolder',figureFolder, ...
    'OutputName','kw_prediction_histories_trackA');


% Save the load-summary numbers used by the two load figures.
kw51_write_csv(fullfile(figureFolder,'kw_load_benchmark_summary.csv'), ...
    loadSummary.headers,loadSummary.rows);

fprintf('\nAll manuscript figures written to:\n%s\n',figureFolder);
end

%% ------------------------------------------------------------------------
function make_fe_mesh(sourceRoot,figureFolder,visible)
% Build the official FE geometry in MATLAB and label the seven plotted
% reference/withheld locations used by the probabilistic model.

stabil=fullfile(sourceRoot,'dependencies','Stabil');
fe=fullfile(sourceRoot,'dependencies','official_fe');
addpath(genpath(stabil));addpath(fe);
old=pwd;cleanup=onCleanup(@()cd(old));cd(fe); %#ok<NASGU>

% Geometry is unchanged by the support choice. plot_elem=0 avoids the
% official routine's own figure; the geometry is redrawn below.
[~,~,~,Nodes,~,~,~,Elements,~]=model_kw51(1,1,0,1,1,0,0,0,[]);

blue=[0.08 0.37 0.51];
gray=[0.72 0.74 0.75];
orange=[0.70 0.38 0.20];

fig=figure('Color','w','Visible',visible,'Units','centimeters', ...
    'Position',[2 2 29.7 15.5]);
layout=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax=nexttile(layout,1);hold(ax,'on');
[nodeID,order]=sort(Nodes(:,1));nodeXYZ=Nodes(order,2:4);

% Shell edges first.
shell=Elements(:,2)==3;
for e=find(shell).'
    ids=Elements(e,5:end);ids=ids(isfinite(ids)&ids>0);
    [ok,idx]=ismember(ids,nodeID);
    idx=idx(ok);
    if numel(idx)<3,continue;end
    xyz=nodeXYZ(idx,:);xyz=[xyz;xyz(1,:)];
    plot3(ax,xyz(:,1),xyz(:,2),xyz(:,3),'-','Color',gray,'LineWidth',0.28);
end
% Beam and truss members.
beam=Elements(:,2)==1 | Elements(:,2)==2;
for e=find(beam).'
    ids=Elements(e,5:6);ids=ids(isfinite(ids)&ids>0);
    [ok,idx]=ismember(ids,nodeID);
    idx=idx(ok);
    if numel(idx)~=2,continue;end
    xyz=nodeXYZ(idx,:);
    plot3(ax,xyz(:,1),xyz(:,2),xyz(:,3),'-','Color',blue,'LineWidth',0.55);
end
axis(ax,'equal');grid(ax,'on');box(ax,'off');view(ax,[-57 19]);
xlabel(ax,'Longitudinal x (m)');ylabel(ax,'y (m)');zlabel(ax,'z (m)');
title(ax,'(a) FE mesh with measurement locations','FontWeight','normal');
set(ax,'FontName','Times New Roman','FontSize',9);

% Measurement markers on the 3-D mesh.
S=sensor_locations();
for k=1:numel(S)
    [mk,mfc,mec]=role_marker(S(k).role,blue,orange);
    plot3(ax,S(k).x,S(k).y,S(k).z,mk,'MarkerSize',6.5, ...
        'MarkerFaceColor',mfc,'MarkerEdgeColor',mec,'LineWidth',1.1);
end

% Girder-axis schematic with readable labels.
ax2=nexttile(layout,2);hold(ax2,'on');
plot(ax2,[0 116.316],[1 1],'k-','LineWidth',1.0);
plot(ax2,[0 116.316],[0 0],'k-','LineWidth',1.0);
text(ax2,-3,1,'A','HorizontalAlignment','right','FontName','Times New Roman');
text(ax2,-3,0,'C','HorizontalAlignment','right','FontName','Times New Roman');
for k=1:numel(S)
    yy=1-(S(k).y/11.75); % A at 1, C at 0
    [mk,mfc,mec]=role_marker(S(k).role,blue,orange);
    plot(ax2,S(k).x,yy,mk,'MarkerSize',7,'MarkerFaceColor',mfc, ...
        'MarkerEdgeColor',mec,'LineWidth',1.2);
    dy=0.18;if yy<0.5,dy=-0.18;end
    text(ax2,S(k).x,yy+dy,sprintf('%s: %s',S(k).key,S(k).channel), ...
        'HorizontalAlignment','center','FontName','Times New Roman','FontSize',8.5, ...
        'Interpreter','none');
end
xlim(ax2,[-4 120]);ylim(ax2,[-0.45 1.45]);yticks(ax2,[]);
xlabel(ax2,'Longitudinal FE coordinate x (m)');
title(ax2,'(b) Girder-axis locations and sensor identities','FontWeight','normal');
box(ax2,'off');set(ax2,'FontName','Times New Roman','FontSize',9);

% Legend dummies.
h1=plot(ax2,NaN,NaN,'s','MarkerFaceColor',blue,'MarkerEdgeColor',blue,'DisplayName','Reference strain');
h2=plot(ax2,NaN,NaN,'o','MarkerFaceColor','w','MarkerEdgeColor',orange,'LineWidth',1.2,'DisplayName','Withheld strain');
h3=plot(ax2,NaN,NaN,'d','MarkerFaceColor','w','MarkerEdgeColor',orange,'LineWidth',1.2,'DisplayName','Withheld displacement');
legend(ax2,[h1 h2 h3],'Location','southoutside','Orientation','horizontal', ...
    'FontName','Times New Roman','FontSize',8.5);

save_figure(fig,fullfile(figureFolder,'kw_fe_mesh'));close_if_hidden(fig,visible);
end

%% ------------------------------------------------------------------------
function make_inputs_figure(runFolder,inputs,figureFolder,visible)
eventID='traindata_20191217_170709';
idx=find(strcmp({inputs.event_id},eventID),1);
assert(~isempty(idx),'Manuscript input event %s is absent from this run.',eventID);
E=load(inputs(idx).input_file);
D=load(fullfile(runFolder,'draws',[eventID '_ensemble_draws.mat']));

blue=[0.08 0.37 0.51];orange=[0.70 0.38 0.20];
fig=figure('Color','w','Visible',visible,'Units','centimeters', ...
    'Position',[2 2 29.7 18.0]);
layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax=nexttile(layout,1);hold(ax,'on');
plot(ax,E.time_s,E.y_ref,'LineWidth',1.0);
xlabel(ax,'Record time (s)');ylabel(ax,'Strain (\muepsilon)');
title(ax,'(a) Reference strain histories','FontWeight','normal');
legend(ax,as_cellstr(E.reference_names),'Interpreter','none','Location','best','FontSize',8);
format_axis(ax);

ax=nexttile(layout,2);hold(ax,'on');
tau=E.rail_offsets_s(:);
for k=1:numel(tau)
    plot(ax,[tau(k) tau(k)],[0 1],'-','Color',blue,'LineWidth',0.85);
end
xlim(ax,[min(0,min(tau)-0.2) max(tau)+0.2]);ylim(ax,[0 1.05]);
xlabel(ax,'Time from first rail peak (s)');ylabel(ax,'Peak indicator');
title(ax,sprintf('(b) Rail-derived candidate timings (%d peaks)',numel(tau)), ...
    'FontWeight','normal');format_axis(ax);

ax=nexttile(layout,3);hold(ax,'on');
loads=double(D.loads_100kN)*100; % kN
[med,lo,hi]=column_quantiles(loads,[0.50 0.05 0.95]);
x=1:size(loads,2);
fill(ax,[x fliplr(x)],[lo fliplr(hi)],[0.82 0.88 0.91], ...
    'EdgeColor','none','DisplayName','5th--95th percentiles');
plot(ax,x,med,'-o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',3.0, ...
    'LineWidth',1.0,'DisplayName','Median');
xlabel(ax,'Candidate axle number');ylabel(ax,'Load (kN)');
title(ax,sprintf('(c) Prior load ensemble: %d candidates',size(loads,1)), ...
    'FontWeight','normal');legend(ax,'Location','best','FontSize',8);format_axis(ax);

ax=nexttile(layout,4);hold(ax,'on');
for k=1:numel(inputs)
    Eall=load(inputs(k).input_file,'event_id','track','rail_offsets_s');
    P=load(fullfile(runFolder,'predictions',[Eall.event_id '_prediction.mat']),'kinematics');
    speed=P.kinematics(1)*3.6;peaks=numel(Eall.rail_offsets_s);
    if strcmp(Eall.track,'A'),col=blue;else,col=orange;end
    plot(ax,speed,peaks,'o','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
end
xlabel(ax,'Reference-inferred entry speed (km/h)');ylabel(ax,'Detected rail peaks');
title(ax,sprintf('(d) All %d accepted passages',numel(inputs)),'FontWeight','normal');
hA=plot(ax,NaN,NaN,'o','Color',blue,'MarkerFaceColor',blue,'DisplayName','Track A');
hB=plot(ax,NaN,NaN,'o','Color',orange,'MarkerFaceColor',orange,'DisplayName','Track B');
legend(ax,[hA hB],'Location','best','FontSize',8);format_axis(ax);

save_figure(fig,fullfile(figureFolder,'kw_inputs'));close_if_hidden(fig,visible);
end

%% ------------------------------------------------------------------------
function S=make_load_benchmark_figures(runFolder,bankFile,inputs,figureFolder,visible)
% Reconstruct exact posterior weights from frozen candidate/paired draws for
% every 12-peak passage. The fixed comparison vector is the measured empty
% Desiro workshop benchmark; it is not an event-specific WIM measurement.

c=load(fullfile(runFolder,'settings.mat'));
B=load(bankFile);
benchmarkAxle=[132.0 131.4 116.8 116.4 94.7 94.6 94.6 94.7 116.4 116.8 131.4 132.0];
benchmarkBogie=sum(reshape(benchmarkAxle,2,6),1);
benchmarkTotal=sum(benchmarkAxle);

names={inputs.event_id};[names,ord]=sort(names);inputs=inputs(ord);
rows={};example=[];
for k=1:numel(inputs)
    E=load(inputs(k).input_file);
    if numel(E.rail_offsets_s)~=12,continue;end
    event=E.event_id;
    D=load(fullfile(runFolder,'draws',[event '_ensemble_draws.mat']));
    P=load(fullfile(runFolder,'predictions',[event '_prediction.mat']), ...
        'kinematics','posterior_ess','track');
    F=struct('kinematics',P.kinematics);
    [nom,delta]=kwpe_simulate(B,E,F,c,D);
    Q=kwpe_condition(nom,delta,E.y_ref,1:3,c,E.time_s);
    assert(abs(Q.posterior_ess-P.posterior_ess)<1e-7, ...
        'Posterior replay mismatch for %s.',event);

    w=Q.weights(:);axle=double(D.loads_100kN)*100;
    [am,al,ah]=weighted_summary(axle,w);
    bogie=reshape(axle,size(axle,1),2,6);bogie=squeeze(sum(bogie,2));
    [bm,bl,bh]=weighted_summary(bogie,w);
    total=sum(axle,2);[tm,tl,th]=weighted_summary(total,w);
    bogieMAE=mean(abs(bm-benchmarkBogie));

    rows(end+1,:)={event,P.track,Q.posterior_ess,tm,tl,th,benchmarkTotal, ...
        mean(abs(am-benchmarkAxle)),bogieMAE}; %#ok<AGROW>
    if isempty(example)
        example=struct('event',event,'track',P.track,'am',am,'al',al,'ah',ah, ...
            'bm',bm,'bl',bl,'bh',bh);
    end
end
assert(~isempty(example),'No accepted 12-peak passages found.');

blue=[0.08 0.37 0.51];orange=[0.70 0.38 0.20];
fig=figure('Color','w','Visible',visible,'Units','centimeters', ...
    'Position',[2 2 29.7 12.8]);
layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(layout,1);hold(ax,'on');
errorbar(ax,1:12,example.am,example.am-example.al,example.ah-example.am,'o', ...
    'Color',blue,'MarkerFaceColor',blue,'LineWidth',1.05,'CapSize',5, ...
    'DisplayName','Posterior mean and 95% interval');
plot(ax,1:12,benchmarkAxle,'s-','Color',orange,'MarkerFaceColor',orange, ...
    'LineWidth',1.05,'DisplayName','Measured empty-train benchmark');
xticks(ax,1:12);xlabel(ax,'Axle number');ylabel(ax,'Load (kN)');
title(ax,'(a) Axle-load comparison','FontWeight','normal');format_axis(ax);

ax=nexttile(layout,2);hold(ax,'on');
errorbar(ax,1:6,example.bm,example.bm-example.bl,example.bh-example.bm,'o', ...
    'Color',blue,'MarkerFaceColor',blue,'LineWidth',1.05,'CapSize',5, ...
    'DisplayName','Posterior mean and 95% interval');
plot(ax,1:6,benchmarkBogie,'s-','Color',orange,'MarkerFaceColor',orange, ...
    'LineWidth',1.05,'DisplayName','Measured empty-train benchmark');
xticks(ax,1:6);xlabel(ax,'Bogie number');ylabel(ax,'Load (kN)');
title(ax,'(b) Bogie-load comparison','FontWeight','normal');format_axis(ax);
legend(ax,'Location','northoutside','FontSize',8);
title(layout,sprintf('%s | Track %s | 12 detected rail peaks', ...
    strrep(example.event,'_','\_'),example.track), ...
    'FontName','Times New Roman','FontSize',10.5,'FontWeight','normal');
save_figure(fig,fullfile(figureFolder,'kw_posterior_load_benchmark'));close_if_hidden(fig,visible);

% All 12-peak passages summary.
eventNames=rows(:,1);n=numel(eventNames);
totalMean=cell2mat(rows(:,4));totalLo=cell2mat(rows(:,5));totalHi=cell2mat(rows(:,6));
bogieMAE=cell2mat(rows(:,9));
fig=figure('Color','w','Visible',visible,'Units','centimeters', ...
    'Position',[2 2 29.7 13.0]);
layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(layout,1);hold(ax,'on');
x=1:n;
errorbar(ax,x,totalMean,totalMean-totalLo,totalHi-totalMean,'o','Color',blue, ...
    'MarkerFaceColor',blue,'MarkerSize',4,'LineWidth',0.9,'CapSize',3);
plot(ax,[0.5 n+0.5],[benchmarkTotal benchmarkTotal],'--','Color',orange,'LineWidth',1.2);
xlim(ax,[0.5 n+0.5]);xticks(ax,x);xticklabels(ax,short_event_labels(eventNames));xtickangle(ax,55);
ylabel(ax,'Total train load (kN)');xlabel(ax,'Passage');
title(ax,'(a) Posterior total load','FontWeight','normal');format_axis(ax);

ax=nexttile(layout,2);bar(ax,x,bogieMAE,'FaceColor',blue,'EdgeColor','none');
xlim(ax,[0.5 n+0.5]);xticks(ax,x);xticklabels(ax,short_event_labels(eventNames));xtickangle(ax,55);
ylabel(ax,'Mean absolute bogie-load difference (kN)');xlabel(ax,'Passage');
title(ax,'(b) Bogie-load difference from benchmark','FontWeight','normal');format_axis(ax);
save_figure(fig,fullfile(figureFolder,'kw_posterior_load_benchmark_17passage_summary'));close_if_hidden(fig,visible);

S.headers={'event_id','track','posterior_ESS','posterior_total_mean_kN', ...
    'posterior_total_q025_kN','posterior_total_q975_kN','benchmark_total_kN', ...
    'axle_mean_MAE_kN','bogie_mean_MAE_kN'};
S.rows=rows;
end

%% ------------------------------------------------------------------------
function make_validation_summary(summary,c,figureFolder,visible)
targets=as_cellstr(c.protected_names);
split=reporting_split(summary);
values=nan(4,3);
for j=1:4
    ix=strcmp(summary(:,1),split)&strcmp(summary(:,2),targets{j})& ...
        abs(cell2mat(summary(:,3))-c.translation_sd_mm)<1e-12;
    k=find(ix,1);
    if ~isempty(k),values(j,:)=[summary{k,5},100*summary{k,8},summary{k,9}];end
end
blue=[0.08 0.37 0.51];
fig=figure('Color','w','Visible',visible,'Units','centimeters', ...
    'Position',[2 2 29.7 11.5]);
layout=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
titles={'(a) Median prediction error','(b) Nominal interval coverage','(c) Interval width'};
ylabels={'Range-normalised RMSE','Mean coverage (%)','Median width / measured range'};
for p=1:3
    ax=nexttile(layout,p);bar(ax,values(:,p),'FaceColor',blue,'EdgeColor','none');hold(ax,'on');
    if p==1
        plot(ax,[0.5 4.5],[0.15 0.15],'k:','LineWidth',1.1);
    elseif p==2
        plot(ax,[0.5 4.5],[95 95],'k--','LineWidth',1.0);ylim(ax,[0 110]);
    end
    xticks(ax,1:4);xticklabels(ax,targets);xtickangle(ax,25);
    ylabel(ax,ylabels{p});title(ax,titles{p},'FontWeight','normal');format_axis(ax);
end
save_figure(fig,fullfile(figureFolder,'kw_ensemble_validation_summary'));close_if_hidden(fig,visible);
end

%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
function root=package_root()
here=fileparts(mfilename('fullpath'));
if exist(fullfile(here,'config.json'),'file')==2
    root=here;return;
end
parent=fileparts(here);
if exist(fullfile(parent,'config.json'),'file')==2
    root=parent;return;
end
error('Could not locate the KW51 package root. Put this function in the package root or functions folder.');
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

function split=reporting_split(summary)
if any(strcmp(summary(:,1),'evaluation'))
    split='evaluation';
elseif any(strcmp(summary(:,1),'development'))
    split='development';
else
    error('No recognised saved reporting split found.');
end
end

function S=sensor_locations()
S=struct( ...
    'key',{'R1','R2','R3','T1','T2','D1','D2'}, ...
    'channel',{'sgBD1011A','sgBD1415A','sgBD2728A','sgBD1718C','sgBD2324C','dABLA','dABLC'}, ...
    'role',{'reference','reference','reference','target_strain','target_strain','target_disp','target_disp'}, ...
    'x',{32.5915,48.3245,99.3080,60.4245,83.7245,115.2570,115.2570}, ...
    'y',{0,0,0,11.75,11.75,0,11.75}, ...
    'z',{0.2719,0.2719,0.2719,0.2719,0.2719,-0.15,-0.15});
end

function [mk,mfc,mec]=role_marker(role,blue,orange)
switch role
    case 'reference'
        mk='s';mfc=blue;mec=blue;
    case 'target_strain'
        mk='o';mfc='w';mec=orange;
    otherwise
        mk='d';mfc='w';mec=orange;
end
end

function [mu,lo,hi]=weighted_summary(X,w)
% Posterior weighted mean and 2.5/97.5 percentiles.
%
% Do not use interp1 directly on cumsum(w): posterior likelihoods can
% underflow to exactly zero (or become too small to change the cumulative
% sum in floating-point arithmetic), which creates repeated CDF abscissae
% and causes MATLAB R2025a interp1 to fail with "Sample points must be
% unique."  Remove non-positive weights and explicitly uniquify the CDF.

X=double(X);
w=double(w(:));

validW=isfinite(w) & w>0;
assert(any(validW),'Posterior weights contain no positive finite values.');

X=X(validW,:);
w=w(validW);
w=w/sum(w);

mu=w.'*X;
lo=zeros(1,size(X,2));
hi=lo;

for j=1:size(X,2)
    validX=isfinite(X(:,j));
    assert(any(validX),'Load candidates contain no finite values.');

    v=X(validX,j);
    wj=w(validX);
    wj=wj/sum(wj);

    [v,idx]=sort(v);
    wj=wj(idx);
    cw=cumsum(wj);
    cw(end)=1; % guard against round-off

    % cumsum can still contain duplicates when a positive weight is below
    % machine resolution relative to the accumulated probability.  Since
    % cw is monotone, retain the last candidate in each repeated-CDF block.
    keep=[diff(cw)>0;true];
    cwu=cw(keep);
    vu=v(keep);

    % Numerical underflow can leave an initial zero-CDF block even though
    % the original weights were positive. Remove it before prepending the
    % explicit zero-probability anchor below.
    positiveCdf=cwu>0;
    cwu=cwu(positiveCdf);
    vu=vu(positiveCdf);
    assert(~isempty(cwu),'Posterior CDF contains no positive probability mass.');

    % Add a zero-probability anchor at the smallest retained candidate.
    xp=[0;cwu(:)];
    vp=[vu(1);vu(:)];

    if numel(xp)==2
        lo(j)=vu(1);
        hi(j)=vu(1);
    else
        lo(j)=interp1(xp,vp,0.025,'linear','extrap');
        hi(j)=interp1(xp,vp,0.975,'linear','extrap');
    end
end
end

function [med,lo,hi]=column_quantiles(X,p)
X=sort(double(X),1);n=size(X,1);
q=zeros(numel(p),size(X,2));
for i=1:numel(p)
    pos=1+(n-1)*p(i);a=floor(pos);b=ceil(pos);f=pos-a;
    if a==b,q(i,:)=X(a,:);else,q(i,:)=(1-f)*X(a,:)+f*X(b,:);end
end
med=q(1,:);lo=q(2,:);hi=q(3,:);
end

function labels=short_event_labels(events)
labels=cell(size(events));
for k=1:numel(events)
    s=char(events{k});
    tok=regexp(s,'traindata_(\d{8})_(\d{6})','tokens','once');
    if isempty(tok),labels{k}=s;continue;end
    d=tok{1};tm=tok{2};
    labels{k}=sprintf('%s/%s %s:%s',d(7:8),d(5:6),tm(1:2),tm(3:4));
end
end

function c=as_cellstr(x)
if iscell(x),c=x(:);elseif isstring(x),c=cellstr(x(:));else,c=cellstr(x);end
end

function format_axis(ax)
grid(ax,'on');box(ax,'off');
set(ax,'FontName','Times New Roman','FontSize',9.5,'LineWidth',0.75, ...
    'TickDir','out');
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

% The FE mesh contains thousands of primitives.  MATLAB R2025a warns
% that vectorizing this particular figure can be very slow and may produce
% unexpected output, so export only the mesh PDF as a high-resolution image.
% All other manuscript plots remain vector PDFs.
useImagePdf=contains(lower(base),'kw_fe_mesh');

if useImagePdf
    exportgraphics(fig,pdfFile, ...
        'ContentType','image', ...
        'Resolution',300, ...
        'BackgroundColor','white');
else
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
end
function close_if_hidden(fig,visible)
if strcmpi(visible,'off'),close(fig);end
end
