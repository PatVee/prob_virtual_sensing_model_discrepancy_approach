function peaks=kw51_peaks(y,height,distance)
% Local maxima and height-priority minimum-distance suppression.
% Plateau maxima use their middle sample, rounded down.
y=y(:);c=[];i=2;
while i<numel(y)
 if y(i)>y(i-1)
  last=i;while last<numel(y)&&y(last+1)==y(i),last=last+1;end
  if last<numel(y)&&y(last)>y(last+1)&&y(i)>=height,c(end+1)=floor((i+last)/2);end %#ok<AGROW>
  i=last+1;
 else,i=i+1;end
end
[~,order]=sort(y(c),'descend');selected=[];
for k=order(:).',if isempty(selected)||all(abs(c(k)-selected)>=distance),selected(end+1)=c(k);end;end %#ok<AGROW>
peaks=sort(selected(:));
end
