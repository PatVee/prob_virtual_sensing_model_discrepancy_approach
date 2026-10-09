function [P,draws]=kwpe_predict(B,E,F,c)
kw51_check_input(E);seed=kwpe_seed(E.event_id,c.seed);rng(seed,'twister');
[nom,delta,draws]=kwpe_simulate(B,E,F,c);
Q=kwpe_condition(nom,delta,E.y_ref,1:3,c,E.time_s);
P=rmfield(Q,'trajectories');
P.raw_target_draws=single(Q.trajectories(:,:,4:7));
P.time_s=E.time_s(:);P.event_id=E.event_id;P.track=E.track;P.labels=B.labels;
P.reference_names=E.reference_names;P.kinematics=F.kinematics;
P.rail_offsets_s=E.rail_offsets_s;P.seed=seed;
P.translation_sd_mm=c.translation_sd_mm;P.translation_sensitivity_mm=c.translation_sensitivity_mm;
P.interval_z=c.interval_z;
P.posterior_load_mean_N=Q.weights.'*draws.loads_100kN*1e5;
P.posterior_stiffness_mean=Q.weights.'*draws.stiffness_states;
P.active_mask=kw51_active_mask(E,F.kinematics);active=P.active_mask;
P.reference_nrmse=sqrt(mean((P.mean(active,1:3)-E.y_ref(active,:)).^2,1))./ ...
    (max(E.y_ref(active,:),[],1)-min(E.y_ref(active,:),[],1));
P.method='KW51 blueprint-BC paired-FE discrepancy covariance with reference-only ensemble conditioning';
P.translation_note='Common-translation prior disabled: blueprint BC anchors longitudinal rigid motion at Brussels.';
if isfield(B,'boundary_description'),P.boundary_description=B.boundary_description;end
draws.seed=seed;draws.event_id=E.event_id;
end
