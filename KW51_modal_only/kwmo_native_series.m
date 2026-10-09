function A=kwmo_native_series(V,j)
% Normalise vectors without flattening a malformed matrix or inventing data.
A=struct('available',false,'time',[],'truth',[],'prediction',[],'lower',[],'upper',[],'mask',[]);
t=V.native_time{j};y=V.native_truth{j};p=V.native_prediction{j};
if isempty(t)||isempty(y)||isempty(p),return;end
assert(isvector(t)&&isvector(y)&&isvector(p)&&numel(t)==numel(y)&&numel(t)==numel(p), ...
 'Unexpected native validation dimensions for target %d.',j);
A.available=true;A.time=t(:);A.truth=y(:);A.prediction=p(:);
lo=V.native_lower95{j,1};hi=V.native_upper95{j,1};
assert(isempty(lo)==isempty(hi),'Only one interval bound was stored.');
if ~isempty(lo)
 assert(isvector(lo)&&isvector(hi)&&numel(lo)==numel(t)&&numel(hi)==numel(t),'Interval dimensions differ from native times.');
 A.lower=lo(:);A.upper=hi(:);
end
if isfield(V,'native_score_mask')&&~isempty(V.native_score_mask{j})
 m=V.native_score_mask{j};assert(isvector(m)&&numel(m)==numel(t));A.mask=logical(m(:));
else,A.mask=false(size(A.time));end
end
