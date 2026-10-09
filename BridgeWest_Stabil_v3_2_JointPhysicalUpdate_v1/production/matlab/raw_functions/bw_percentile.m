function y=bw_percentile(x,p)
% NumPy default linear percentile, avoiding a Statistics Toolbox dependency.
x=sort(x(:)); a=1+(numel(x)-1)*p/100; lo=floor(a); hi=ceil(a);
y=x(lo)+(a-lo)*(x(hi)-x(lo));
end
