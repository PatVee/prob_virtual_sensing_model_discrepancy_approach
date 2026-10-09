function [P,draws,D]=kwmo_predict(B,E,F,c,S,checkConvergence)
kw51_check_input(E);seed=kwpe_seed(E.event_id,c.seed);
rng(seed,'twister');loads=kwpe_load_draws(c.candidate_count,numel(E.rail_offsets_s),c);
pairloads=kwpe_load_draws(c.paired_count,numel(E.rail_offsets_s),c);
nom=kwmo_nominal(B,E,F,loads);
rng(seed+301,'twister');[delta,modal]=kwmo_modal_delta(B,E,F,c,pairloads,S,'PointsPerCycle',c.points_per_cycle);
convergence=struct('checked',false,'relative_error_by_output',nan(1,7));
if checkConvergence
 rng(seed+301,'twister');fine=kwmo_modal_delta(B,E,F,c,pairloads,S,'PointsPerCycle',2*c.points_per_cycle);
 for o=1:7
  convergence.relative_error_by_output(o)=norm(delta(:,:,o)-fine(:,:,o),'fro')/max(norm(fine(:,:,o),'fro'),1e-12);
 end
 convergence.checked=true;
 assert(all(convergence.relative_error_by_output<c.integration_relative_tolerance), ...
  'Modal integration convergence failed; increase PointsPerCycle.');
end
rng(seed+907,'twister');Q=kwpe_condition(nom,delta,E.y_ref,1:3,c,E.time_s);
D=kwmo_decompose(nom,delta,E.y_ref,Q,c);D.time_s=E.time_s(:);D.event_id=E.event_id;D.labels=B.labels;
P=rmfield(Q,'trajectories');P.raw_target_draws=single(Q.trajectories(:,:,4:7));
P.time_s=E.time_s(:);P.event_id=E.event_id;P.track=E.track;P.labels=B.labels;
P.reference_names=E.reference_names;P.kinematics=F.kinematics;P.rail_offsets_s=E.rail_offsets_s;P.seed=seed;
P.translation_sd_mm=0;P.translation_sensitivity_mm=0;P.interval_z=c.interval_z;
P.posterior_load_mean_N=Q.weights.'*loads*1e5;P.baseline_stiffness_scale=1;
P.active_mask=kw51_active_mask(E,F.kinematics);a=P.active_mask;
P.reference_nrmse=sqrt(mean((P.mean(a,1:3)-E.y_ref(a,:)).^2,1))./ ...
 (max(E.y_ref(a,:),[],1)-min(E.y_ref(a,:),[],1));
P.method=c.method;P.boundary_description=B.boundary_description;P.integration_convergence=convergence;
P.modal_parameter_note=S.policy;P.nominal_model='Fixed blueprint quasi-static response';
draws=struct('loads_100kN',loads,'paired_loads_100kN',pairloads,'modal_draws',modal, ...
 'seed',seed,'event_id',E.event_id,'baseline_stiffness_scale',1);
end
