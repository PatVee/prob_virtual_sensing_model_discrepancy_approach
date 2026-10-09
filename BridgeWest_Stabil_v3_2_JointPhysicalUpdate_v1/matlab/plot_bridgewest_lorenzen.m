function figs=plot_bridgewest_lorenzen(fe,outputDir)
% Plotting is separate from the solver. Titles state the provisional status.
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1||isempty(fe)
    z=load(fullfile(root,'results','matlab','bridgewest_lorenzen_fe.mat'),'fe');fe=z.fe;
end
if nargin<2,outputDir=fullfile(root,'figures','matlab');end
if ~isfolder(outputDir),mkdir(outputDir);end
figs.geometry=figure('Color','w','Position',[100 100 1200 650]);ax=axes;hold(ax,'on');
patch(ax,'Vertices',fe.nodes,'Faces',fe.shells(:,1:4),'FaceColor',[.70 .78 .83],'EdgeColor',[.35 .4 .45],'FaceAlpha',.75);
for i=1:size(fe.beams,1),p=fe.nodes(fe.beams(i,1:2),:);plot3(ax,p(:,1),p(:,2),p(:,3),'k-');end
support=unique(floor((fe.fixed-1)/6)+1);p=fe.nodes(support,:);scatter3(ax,p(:,1),p(:,2),p(:,3),60,'r','filled','^');
q=find(fe.foundation>0);p=fe.nodes(unique(floor((q-1)/6)+1),:);scatter3(ax,p(:,1),p(:,2),p(:,3),18,[0 .5 .2],'filled');
axis(ax,'equal');view(ax,[-25 20]);xlabel(ax,'Model x (m)');ylabel(ax,'y (m)');zlabel(ax,'z (m)');
title(ax,'Lorenzen-inspired Stabil reconstruction — provisional section dimensions');grid(ax,'on');
figs.modes=figure('Color','w','Position',[100 100 1400 900]);tl=tiledlayout(2,3,'TileSpacing','compact');
for i=1:5
    ax=nexttile(tl);hold(ax,'on');j=fe.comparison.fe_index(i);
    if isnan(j)
        axis(ax,'off');text(ax,.05,.7,sprintf('Thesis mode %d: %.2f Hz\n%s\nUnassigned',i,fe.comparison.frequency_Hz(i),fe.comparison.type(i)),'Units','normalized');continue
    end
    U=reshape(fe.P(:,j),6,[])';fac=.35/max(abs(U(:,1:3)),[],'all');
    patch(ax,'Vertices',fe.nodes+fac*U(:,1:3),'Faces',fe.shells(:,1:4),'FaceVertexCData',U(:,3), ...
        'FaceColor','interp','EdgeColor','none');axis(ax,'equal');view(ax,[-25 20]);grid(ax,'on');
    title(ax,sprintf('Thesis %d: %s\nFE %.2f Hz | observed %.2f Hz',i,fe.comparison.type(i),fe.freq(j),fe.comparison.frequency_Hz(i)));
    xlabel(ax,'x (m)');ylabel(ax,'y (m)');zlabel(ax,'z (m)');
end
ax=nexttile(tl);axis(ax,'off');text(ax,.03,.9,{'Family-based candidates only.','No numerical ambient mode vectors supplied.','Mode 5 is undefined in the thesis.','Geometry has not passed field validation.'},'Units','normalized','VerticalAlignment','top');
figs.profiles=figure('Color','w','Position',[100 100 1200 750]);tl=tiledlayout(2,2,'TileSpacing','compact');
for i=1:4
    ax=nexttile(tl);j=fe.comparison.fe_index(i);
    if isnan(j),axis(ax,'off');title(ax,'No qualifying candidate');continue;end
    A=fe.P(6*(fe.top(1,:)-1)+3,j);B=fe.P(6*(fe.top(2,:)-1)+3,j);scale=max(abs([A;B]));
    plot(ax,fe.x,A/scale,'k-',fe.x,B/scale,'--','LineWidth',1.6);grid(ax,'on');xlabel(ax,'Model x (m)');ylabel(ax,'Normalised vertical mode');
    title(ax,sprintf('%s: %.2f Hz',fe.comparison.type(i),fe.freq(j)));legend(ax,{'BA','BB'},'Location','best');
end
names=fieldnames(figs);
for i=1:numel(names)
    exportgraphics(figs.(names{i}),fullfile(outputDir,['bw_lorenzen_' names{i} '.png']),'Resolution',220);
    exportgraphics(figs.(names{i}),fullfile(outputDir,['bw_lorenzen_' names{i} '.pdf']),'ContentType','vector');
end
end
