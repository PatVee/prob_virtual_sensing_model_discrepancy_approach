function plot_kw51_modal_only(out,varargin)
% Plot measured histories, modal-only posterior mean and nominal 95% interval.
ip=inputParser;addParameter(ip,'Events','all');addParameter(ip,'Visible','off');parse(ip,varargin{:});
if ischar(out)||isstring(out),L=load(fullfile(char(out),'run_output.mat'));out=L.out;end
kwmo_setup(out.original_root);kwpe_verify_freeze(out.run_dir);
I=load(fullfile(out.run_dir,'input_catalog.mat'));dest=fullfile(out.run_dir,'figures');if ~isfolder(dest),mkdir(dest);end
requested=ip.Results.Events;if ischar(requested),requested={requested};elseif isstring(requested),requested=cellstr(requested);end
blue=[.1 .35 .7];
for k=1:numel(I.inputs)
 id=I.inputs(k).event_id;if ~any(strcmp(requested,'all'))&&~any(strcmp(requested,id)),continue;end
 vf=fullfile(out.run_dir,'validation',[id '_validation.mat']);assert(isfile(vf),'Run target validation first.');V=load(vf);
 fig=figure('Color','w','Visible',ip.Results.Visible,'Units','centimeters','Position',[2 2 29 17]);
 tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
 for j=1:4
  ax=nexttile(tl,j);hold(ax,'on');title(ax,V.labels{j},'Interpreter','none','FontWeight','normal');
  A=kwmo_native_series(V,j);
  if ~A.available
   text(ax,.5,.5,{'Validation unavailable','See validation_skipped.csv'},'Units','normalized','HorizontalAlignment','center');
   axis(ax,'off');continue;
  end
  h=gobjects(0);names={};
  if ~isempty(A.lower)
   valid=isfinite(A.time)&isfinite(A.lower)&isfinite(A.upper);
   if any(valid)
    t=A.time(valid);lo=A.lower(valid);hi=A.upper(valid);
    h(end+1)=fill(ax,[t;flipud(t)],[lo;flipud(hi)],blue,'FaceAlpha',.18,'EdgeColor','none');names{end+1}='Nominal 95% interval'; %#ok<AGROW>
   end
  end
  h(end+1)=plot(ax,A.time,A.truth,'k','LineWidth',1.1);names{end+1}='Measured'; %#ok<AGROW>
  h(end+1)=plot(ax,A.time,A.prediction,'Color',blue,'LineWidth',1.2);names{end+1}='Predicted mean'; %#ok<AGROW>
  xlabel(ax,'Time (s)');if j<=2,ylabel(ax,'Strain (\mu\epsilon)');else,ylabel(ax,'Longitudinal displacement (mm)');end
  grid(ax,'on');box(ax,'off');set(ax,'FontSize',12);
  legend(ax,h,names,'Location','best','Box','off');
 end
 title(tl,id,'Interpreter','none');
 exportgraphics(fig,fullfile(dest,[id '_modal.pdf']),'ContentType','vector');
 exportgraphics(fig,fullfile(dest,[id '_modal.png']),'Resolution',200);
 if strcmp(ip.Results.Visible,'off'),close(fig);end
end
end
