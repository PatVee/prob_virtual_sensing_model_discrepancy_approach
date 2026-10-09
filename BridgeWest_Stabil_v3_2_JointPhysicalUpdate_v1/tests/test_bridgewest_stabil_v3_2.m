function test_bridgewest_stabil_v3_2(outputRoot)
%TEST_BRIDGEWEST_STABIL_V3_2 Native-MATLAB selection/data tests; no plots needed.
% Call after adding production/matlab to the MATLAB path.
a=bw_select_saved_events(outputRoot,'Event','all');
b=bw_select_saved_events(outputRoot,'Campaign',2022);
c=bw_select_saved_events(outputRoot,'Campaign',2025);
assert(numel(a)==36&&numel(b)==12&&numel(c)==24);
assert(isequal(bw_select_saved_events(outputRoot,'Event','2025_run22'),{'2025_run22'}));
assert(isequal(bw_select_saved_events(outputRoot,'Event','2022_run20'),{'2022_run20'}));
x=bw_select_saved_events(outputRoot,'Campaign',2025,'Runs',[20 22 40]);
assert(isequal(x,{'2025_run20';'2025_run22';'2025_run40'}));
x=bw_select_saved_events(outputRoot,'Event',{'2022_run20','2025_run22'});assert(numel(x)==2);
rejected=false;try,bw_select_saved_events(outputRoot,'Event','2022_run22');catch ME,rejected=strcmp(ME.identifier,'BridgeWest:MissingEvent');end
assert(rejected,'The configured 2022 set has no run22; it must not be silently substituted.');
T=readtable(fullfile(outputRoot,'tables','target_assessment_all.csv'),'TextType','string');n=0;
for i=1:numel(a)
    D=bw_stabil_strain_plot_data(outputRoot,a{i});assert(isequal(D.labels,{'A058','A070','B031','B081'}));
    assert(size(D.mean,2)==4&&size(D.mean,1)==numel(D.time_s));
    for j=1:4
        q=string(T.event)==a{i}&string(T.target)==D.labels{j};assert(nnz(q)==1);
        assert(logical(T.primary_score(q))==D.primary_score(j));
        if D.primary_score(j)
            assert(abs(D.nrmse(j)-T.nrmse(q))<1e-10);
            assert(abs(D.coverage(j)-T.coverage(q))<1e-10);
            assert(abs(D.correlation(j)-T.correlation(q))<1e-10);n=n+1;
        end
    end
end
assert(n==132);W=bw_weighed_bogie_loads(2022);assert(abs(W.total_load_kN-2206.89832265)<1e-8);
W=bw_weighed_bogie_loads(2025);assert(abs(W.total_load_kN-2197.38587215)<1e-8);
fprintf('PASS: 36 event selections, 144 histories, 132 primary scores, two weighing totals.\n');
end
