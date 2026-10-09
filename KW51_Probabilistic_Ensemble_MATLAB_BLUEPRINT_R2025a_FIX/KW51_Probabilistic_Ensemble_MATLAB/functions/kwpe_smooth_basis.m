function psi=kwpe_smooth_basis(t,maxFrequency)
t=t(:);duration=t(end)-t(1);assert(duration>0);
freq=(1:floor(duration*maxFrequency))/duration;
psi=[ones(numel(t),1),cos(2*pi*t*freq),sin(2*pi*t*freq)]/sqrt(1+numel(freq));
end
