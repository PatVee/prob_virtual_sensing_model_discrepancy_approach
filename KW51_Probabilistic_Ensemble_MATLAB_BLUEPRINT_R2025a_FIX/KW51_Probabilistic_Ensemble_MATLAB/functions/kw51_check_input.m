function kw51_check_input(E)
refs={'sgBD1011A';'sgBD1415A';'sgBD2728A'};
allowed={'time_s','y_ref','reference_names','rail_offsets_s','rail_first_time_s', ...
 'rail_relative_load','track','event_id','source_member','strain_start_sdn','raw_archive'};
assert(all(ismember(fieldnames(E),allowed)),'KW51:ForbiddenInput','Undeclared or protected input field.');
assert(isequal(E.reference_names(:),refs),'KW51:ReferenceList','Reference whitelist mismatch.');
assert(size(E.y_ref,2)==3 && all(isfinite(E.y_ref(:))),'Invalid reference data.');
assert(numel(E.time_s)==size(E.y_ref,1) && all(diff(E.time_s)>0),'Invalid time grid.');
assert(all(diff(E.rail_offsets_s)>0) && E.rail_offsets_s(1)==0,'Invalid axle times.');
assert(any(strcmp(E.track,{'A','B'})),'Invalid track.');
end
