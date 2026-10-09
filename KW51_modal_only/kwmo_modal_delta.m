function [delta,D]=kwmo_modal_delta(B,E,F,c,pairedLoads,S,varargin)
% Paired frequency/damping perturbations about the unscaled blueprint FE baseline.
ip=inputParser;addParameter(ip,'PointsPerCycle',40);addParameter(ip,'MinimumPaddingSeconds',5);addParameter(ip,'SampleIndices',[]);parse(ip,varargin{:});
p=ip.Results;assert(p.PointsPerCycle>=20&&p.MinimumPaddingSeconds>=0);
assert(c.baseline_stiffness_scale==1,'The modal-only workflow uses the unscaled FE baseline.');
t=E.time_s(:);kin=F.kinematics(:).';modes=S.fe_modes;nm=numel(modes);
f0=reshape(B.frequency_Hz(modes),1,[]);z0=S.median_damping;
% M=size(pairedLoads,1);selected=randi(size(S.frequency_ratio,1),M,1);
% f=f0.*S.frequency_ratio(selected,:);z=S.damping_ratio(selected,:);
M = size(pairedLoads,1);selected = p.SampleIndices;
if isempty(selected);selected = randi(size(S.frequency_ratio,1),M,1);else;selected = selected(:);end
assert(numel(selected)==M&&all(selected>=1&selected<=size(S.frequency_ratio,1)&selected==floor(selected)),'Invalid modal sample indices.');
f=f0.*S.frequency_ratio(selected,:);z=S.damping_ratio(selected,:);
maxf=max([f(:);f0(:)]);h=min(median(diff(t))/4,1/(p.PointsPerCycle*maxf));
% Start before the leading axle reaches even the FE force-grid boundary.
assert(kin(1)>0,'Positive speed magnitude required.');
pad=max(p.MinimumPaddingSeconds,abs(min(B.xgrid))/kin(1)+2);
ts=min(t(1),kin(2)-pad);n=ceil((t(end)-ts)/h);tt=linspace(ts,t(end),n+1).';h=tt(2)-tt(1);
tau=E.rail_offsets_s(:).';dt=tt-kin(2);vr=kin(1)+kin(4)*(E.rail_first_time_s-kin(2));
x=kin(1)*dt+.5*kin(4)*dt.^2-(vr*tau+.5*kin(4)*tau.^2);
if kin(3)<0,x=116.316-x;end
assert(all(x(1,:)<min(B.xgrid)|x(1,:)>max(B.xgrid)), ...
 'Force already acts at padded initial time; increase padding/check kinematics.');
speed2=(kin(1)+kin(4)*dt).^2;tr=1+strcmp(E.track,'B');
forceBasis=zeros(numel(tt),nm,numel(tau));
for a=1:numel(tau)
 v0=interp1(B.xgrid,B.modal_force0(:,modes,tr),x(:,a),'linear',0);
 v2=interp1(B.xgrid,B.modal_force2(:,modes,tr),x(:,a),'linear',0);
 forceBasis(:,:,a)=(v0+speed2.*v2)*1e5;
end
delta=zeros(M,numel(t),7);
for m=1:M
 force=reshape(reshape(forceBasis,[],numel(tau))*pairedLoads(m,:).',numel(tt),nm);
 q0=kwmo_newmark(force,h,f0,z0);qp=kwmo_newmark(force,h,f(m,:),z(m,:));
 dy=(qp-q0)*B.output_modes(:,modes).';
 delta(m,:,:)=reshape(interp1(tt,dy,t,'linear'),1,numel(t),7);
end
D=struct('sample_indices',selected,'frequency_Hz',f,'damping_ratio',z, ...
 'nominal_frequency_Hz',f0,'nominal_damping_ratio',z0,'fe_modes',modes, ...
 'integration_dt_s',h,'integration_start_s',ts,'points_per_cycle',p.PointsPerCycle, ...
 'paired_loads_100kN',pairedLoads,'modal_source_sha256',S.source_sha256);
end
