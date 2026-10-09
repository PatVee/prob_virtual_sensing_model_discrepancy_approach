function [p,diagout]=bw_predict_continuous_joint(B0,components,yref,refs,c)
%BW_PREDICT_CONTINUOUS_JOINT Bounded load inference with joint discrepancy.
% Only permitted reference strain yref is accepted. Static weighing, PLS and
% protected target strains are not function arguments.
assert(numel(refs)==size(yref,2));nb=8;
pbar=repmat((c.load_lower_kN+c.load_upper_kN)/2,nb,1);
trace=struct('iteration',{},'max_bogie_change_kN',{},'condition_number',{},'posterior_mean_kN',{});
for itn=1:c.fixed_point_max_iterations
    D=bw_joint_discrepancy_weighted(components,pbar.'*1000,refs,c);
    q=bw_continuous_load_posterior(B0{1}(:,refs,:),yref,D,c,'Seed',c.gibbs_seed, ...
        'Burn',c.gibbs_burn_iteration,'Keep',c.gibbs_keep_iteration);
    raw=q.mean;proposed=(1-c.fixed_point_relaxation)*pbar+c.fixed_point_relaxation*raw;
    change=max(abs(proposed-pbar));
    trace(end+1)=struct('iteration',itn,'max_bogie_change_kN',change, ...
        'condition_number',q.condition_number,'posterior_mean_kN',raw.'); %#ok<AGROW>
    pbar=proposed;if change<c.fixed_point_tolerance_kN,break;end
end
for pass=1:2
    D=bw_joint_discrepancy_weighted(components,pbar.'*1000,refs,c);
    q=bw_continuous_load_posterior(B0{1}(:,refs,:),yref,D,c,'Seed',c.gibbs_seed+100+pass-1, ...
        'Burn',c.gibbs_burn_final,'Keep',c.gibbs_keep_final);
    pbar=q.mean;
end
p.load_mean_kN=q.mean;p.load_covariance_kN2=q.covariance;p.load_sd_kN=sqrt(diag(q.covariance));
p.load_correlation=q.covariance./sqrt(max(diag(q.covariance),eps)*max(diag(q.covariance),eps).');
p.load_condition_number=q.condition_number;p.load_time_indices=q.time_indices;p.iteration_trace=trace;
keys={'strain','disp_BA','disp_BB'};BrAll=B0{1}(:,refs,:);
for k=1:3
    T=size(B0{k},1);O=size(B0{k},2);mn=zeros(T,O);sd=zeros(T,O);
    for it=1:T
        B=reshape(B0{k}(it,:,:),O,nb)*1000;Br=reshape(BrAll(it,:,:),numel(refs),nb)*1000;
        G=D.conditional_gain{k}(:,:,it);A=B-G*Br;
        cc=D.mean{k}(it,:).'+G*(yref(it,:).'-D.mean{1}(it,refs).');
        mn(it,:)=(cc+A*q.mean).';
        vload=sum((A*q.covariance).*A,2);sd(it,:)=sqrt(max(vload+D.conditional_variance{k}(it,:).',0)).';
    end
    p.([keys{k} '_mean'])=mn;p.([keys{k} '_sd'])=sd;
    p.([keys{k} '_lo'])=mn-c.interval_z*sd;p.([keys{k} '_hi'])=mn+c.interval_z*sd;
end
diagout=struct('discrepancy',D,'load_posterior',q,'fixed_point_load_kN',pbar,'trace',{trace});
end
