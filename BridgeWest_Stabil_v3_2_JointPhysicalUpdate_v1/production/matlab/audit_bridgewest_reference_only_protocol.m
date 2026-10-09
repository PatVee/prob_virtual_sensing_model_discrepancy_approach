function audit=audit_bridgewest_reference_only_protocol(outputRoot)
%AUDIT_BRIDGEWEST_REFERENCE_ONLY_PROTOCOL Verify the saved 36-passage firewall.
root=fileparts(mfilename('fullpath'));pkg=fileparts(fileparts(root));
c=jsondecode(fileread(fullfile(pkg,'production','config','prediction_config.json')));
ev=jsondecode(fileread(fullfile(pkg,'production','config','event_sets.json')));
refs=cellstr(string(c.measured_reference_channels));targets=cellstr(string(c.measured_target_channels));
expected=[arrayfun(@(r)sprintf('2022_run%02d',r),double(ev.runs_2022(:)).','UniformOutput',false), ...
          arrayfun(@(r)sprintf('2025_run%02d',r),double(ev.runs_2025(:)).','UniformOutput',false)];
assert(numel(expected)==36,'Configured default event set is not 36 passages.');
rows=struct('event',{},'prediction_exists',{},'permitted_input_exists',{},'reference_count',{}, ...
    'reference_names_exact',{},'target_leakage',{},'prediction_declares_no_target_use',{}, ...
    'validation_target_file_exists',{});
for i=1:numel(expected)
    eid=expected{i};pf=fullfile(outputRoot,'predictions',[eid '.mat']);rf=fullfile(outputRoot,'permitted_inputs',[eid '.mat']);
    row=struct('event',eid,'prediction_exists',exist(pf,'file')==2,'permitted_input_exists',exist(rf,'file')==2, ...
        'reference_count',0,'reference_names_exact',false,'target_leakage',false, ...
        'prediction_declares_no_target_use',false,'validation_target_file_exists',exist(fullfile(outputRoot,'validation_targets',[eid '.mat']),'file')==2);
    if row.permitted_input_exists
        R=load(rf);if isfield(R,'reference_names'),rn=cellstr(string(R.reference_names));else,rn={};end
        if isfield(R,'reference_strain_microstrain'),row.reference_count=size(R.reference_strain_microstrain,2);end
        row.reference_names_exact=isequal(rn(:),refs(:));
        row.target_leakage=any(ismember(rn,targets));
    end
    if row.prediction_exists
        P=load(pf,'protected_targets_used_in_inference','target_measurements_loaded_during_inference');
        a=isfield(P,'protected_targets_used_in_inference')&&~logical(P.protected_targets_used_in_inference);
        b=isfield(P,'target_measurements_loaded_during_inference')&&~logical(P.target_measurements_loaded_during_inference);
        row.prediction_declares_no_target_use=a&&b;
    end
    rows(end+1)=row; %#ok<AGROW>
end
audit=struct2table(rows);tdir=fullfile(outputRoot,'tables');if ~exist(tdir,'dir'),mkdir(tdir);end
writetable(audit,fullfile(tdir,'reference_only_protocol_audit.csv'));
complete=audit.prediction_exists&audit.permitted_input_exists&audit.reference_count==2&audit.reference_names_exact&~audit.target_leakage&audit.prediction_declares_no_target_use;
fprintf('Reference-only protocol audit: %d/%d configured passages pass saved-file checks.\n',sum(complete),height(audit));
if all(audit.prediction_exists),assert(all(complete),'One or more saved 36-passage predictions fail the reference-only audit.');end
end
