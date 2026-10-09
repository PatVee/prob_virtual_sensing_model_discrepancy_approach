function loads=kwpe_load_draws(count,axles,c)
lo=c.load_bounds_kN(1)/100;hi=c.load_bounds_kN(2)/100;
common=lo+(hi-lo)*rand(count,1);
loads=min(max(common.*exp(c.axle_relative_log_sd*randn(count,axles)),lo),hi);
end
