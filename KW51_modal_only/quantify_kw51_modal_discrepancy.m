function quantify_kw51_modal_discrepancy(out,varargin)
% Inferred correction vs posterior prediction; not identification of true error.
ip=inputParser;addParameter(ip,'Plot',false);addParameter(ip,'Visible','off');parse(ip,varargin{:});
if ischar(out)||isstring(out),L=load(fullfile(char(out),'run_output.mat'));out=L.out;end
kwmo_setup(out.original_root);kwpe_verify_freeze(out.run_dir);
I=load(fullfile(out.run_dir,'input_catalog.mat'));dest=fullfile(out.run_dir,'discrepancy_diagnostics');if ~isfolder(dest),mkdir(dest);end
rows={};
for k=1:numel(I.inputs)
 id=I.inputs(k).event_id;D=load(fullfile(out.run_dir,'decomposition',[id '_decomposition.mat']));
 P=load(fullfile(out.run_dir,'predictions',[id '_prediction.mat']));
 V=load(fullfile(out.run_dir,'validation',[id '_validation.mat']));
 t=D.time_s(:);fs=1/median(diff(t));base=t>=0&t<1;
 if ip.Results.Plot
  fig=figure('Color','w','Visible',ip.Results.Visible,'Units','centimeters','Position',[2 2 29 17]);
  tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
 end
 for j=1:4
  o=j+3;A=kwmo_native_series(V,j);
  if ip.Results.Plot
   ax=nexttile(tl,j);hold(ax,'on');title(ax,V.labels{j},'Interpreter','none','FontWeight','normal');
  end
  if ~A.available
   if ip.Results.Plot,text(ax,.5,.5,'Validation unavailable','Units','normalized','HorizontalAlignment','center');axis(ax,'off');end
   continue;
  end
  cut=V.validation_cutoff_Hz(j);
  fe=kw51_filter(D.fe_mean(:,o),cut,fs);dc=kw51_filter(D.correction_mean(:,o),cut,fs);
  corrected=fe+dc;offset=median(corrected(base));fe=fe-median(fe(base));dc=corrected-offset-fe;
  f=interp1(t,fe,A.time);d=interp1(t,dc,A.time);pred=f+d;
  assert(max(abs(pred-A.prediction))<1e-8*max(1,max(abs(pred))),'Filtered decomposition does not reproduce the validated prediction.');
  mask=A.mask&isfinite(f)&isfinite(d)&isfinite(A.truth);
  rmse0=NaN;rmse1=NaN;corrRMS=NaN;rel=NaN;imp=NaN;
  if any(mask)
   rmse0=sqrt(mean((f(mask)-A.truth(mask)).^2));rmse1=sqrt(mean((pred(mask)-A.truth(mask)).^2));
   corrRMS=sqrt(mean(d(mask).^2));feRMS=sqrt(mean(f(mask).^2));
   if feRMS>1e-12,rel=corrRMS/feRMS;end;if rmse0>1e-12,imp=1-rmse1/rmse0;end
  end
  a=P.active_mask;
  rows(end+1,:)={id,V.labels{j},sum(mask),rmse0,rmse1,imp,corrRMS,rel, ...
   mean(D.fe_variance(a,o)),mean(D.correction_between_variance(a,o)), ...
   mean(D.twice_fe_correction_covariance(a,o)),mean(D.conditional_discrepancy_variance(a,o)), ...
   mean(D.output_noise_variance(a,o)),mean(D.total_variance(a,o))}; %#ok<AGROW>
  kw51_write_csv(fullfile(dest,[id '_' V.labels{j} '.csv']), ...
   {'time_s','measured','posterior_FE','inferred_correction','corrected'},num2cell([A.time,A.truth,f,d,pred]));
  if ip.Results.Plot
   plot(ax,A.time,A.truth,'k','LineWidth',1);plot(ax,A.time,f,'--','Color',[.4 .4 .4],'LineWidth',1);
   plot(ax,A.time,pred,'Color',[.1 .35 .7],'LineWidth',1.2);plot(ax,A.time,d,':','Color',[.85 .35 .12],'LineWidth',1.2);
   xlabel(ax,'Time (s)');if j<=2,ylabel(ax,'Strain (\mu\epsilon)');else,ylabel(ax,'Longitudinal displacement (mm)');end
   box(ax,'off');grid(ax,'on');set(ax,'FontSize',12);
   legend(ax,{'Measured','Posterior-weighted FE','Corrected prediction','Inferred correction'},'Box','off','Location','best');
  end
 end
 if ip.Results.Plot
  title(tl,id,'Interpreter','none');
  exportgraphics(fig,fullfile(dest,[id '_decomposition.pdf']),'ContentType','vector');
  exportgraphics(fig,fullfile(dest,[id '_decomposition.png']),'Resolution',200);
  if strcmp(ip.Results.Visible,'off'),close(fig);end
 end
end
kw51_write_csv(fullfile(dest,'discrepancy_summary.csv'), ...
 {'event','target','n_samples','posterior_FE_RMSE','corrected_RMSE','relative_RMSE_reduction', ...
 'correction_RMS','correction_to_FE_RMS','raw_load_variance','raw_correction_between_variance', ...
 'raw_twice_FE_correction_covariance','raw_conditional_discrepancy_variance','raw_output_noise_variance','raw_total_variance'},rows);
end
