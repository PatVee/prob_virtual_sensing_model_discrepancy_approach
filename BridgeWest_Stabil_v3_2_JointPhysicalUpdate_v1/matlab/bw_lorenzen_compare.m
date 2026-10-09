function T=bw_lorenzen_compare(spectrum,targets)
% Unique family/lobe-based screening; never claim MAC from frequencies alone.
family=["bending";"torsion";"bending";"torsion"];lobes=[1;1;2;2];
id=nan(5,1);freq=id;err=id;status=repmat("unassigned",5,1);used=[];
for i=1:4
    rows=find(spectrum.family==family(i) & spectrum.lobes==lobes(i) & ~ismember(spectrum.fe_index,used));
    if ~isempty(rows)
        [~,k]=min(abs(log(spectrum.frequency_Hz(rows)/targets.frequency_Hz(i))));r=rows(k);
        id(i)=spectrum.fe_index(r);freq(i)=spectrum.frequency_Hz(r);err(i)=100*(freq(i)/targets.frequency_Hz(i)-1);
        status(i)="family-based candidate; measured MAC unavailable";used(end+1)=id(i); %#ok<AGROW>
    end
end
status(5)="undefined in thesis; no automatic assignment";
T=targets;T.fe_index=id;T.fe_frequency_Hz=freq;T.frequency_error_percent=err;T.assignment_status=status;
end
