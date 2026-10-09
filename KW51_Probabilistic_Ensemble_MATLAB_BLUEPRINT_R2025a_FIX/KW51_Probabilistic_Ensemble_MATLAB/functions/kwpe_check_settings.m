function kwpe_check_settings(c)
assert(isscalar(c.seed)&&isfinite(c.seed)&&c.seed>=0&&c.seed==floor(c.seed),'Invalid seed.');
assert(c.candidate_count>=3&&c.candidate_count==floor(c.candidate_count)&& ...
    mod(c.candidate_count,numel(c.stiffness_states))==0,'CandidateCount must be divisible by the number of stiffness states.');
assert(c.paired_count>=4&&c.paired_count==floor(c.paired_count),'PairedCount must be an integer >= 4.');
assert(c.posterior_draw_count>=16&&c.posterior_draw_count==floor(c.posterior_draw_count),'PosteriorDrawCount must be an integer >= 16.');
assert(isscalar(c.translation_sd_mm)&&abs(c.translation_sd_mm)<1e-12, ...
    'Blueprint BC requires translation_sd_mm = 0.');
assert(numel(c.translation_sensitivity_mm)==1&&abs(c.translation_sensitivity_mm(1))<1e-12, ...
    'Blueprint BC requires translation_sensitivity_mm = 0 only.');
assert(all(c.stiffness_states>0)&&c.pair_nominal_stiffness>0,'Stiffness scales must be positive.');
assert(isequal(c.reference_names(:),{'sgBD1011A';'sgBD1415A';'sgBD2728A'}),'Reference whitelist changed.');
assert(isequal(c.protected_names(:),{'sgBD1718C';'sgBD2324C';'dABLA';'dABLC'}),'Protected-output list changed.');
end
