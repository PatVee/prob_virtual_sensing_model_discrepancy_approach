function F=kw51_fit(B,E,refs,forceScale)
% Fresh speed, entry, acceleration and load inference. Reads no results file.
kw51_check_input(E);if nargin<3,refs=1:3;end;if nargin<4,forceScale=8;end
t=E.time_s(1:5:end);t=t(:);y=E.y_ref(1:5:end,refs);
[~,W,~]=kw51_force_map(B,refs,forceScale);tau=E.rail_offsets_s(:);
desiro=[0;2.3;18.6;20.9;26.1;28.4;44.7;47;52.2;54.5;70.8;73.1];matched=false;railseed=[];
if numel(tau)==12
 X=[tau,.5*tau.^2];va=X\desiro;err=sqrt(mean((X*va-desiro).^2));
 if va(1)>8 && va(1)<44.44 && abs(va(2))<.8 && err<.35,matched=true;railseed=[va(1),118,va(2)];end
end
candidates=zeros(0,5);
for direction=[-1,1]
 coarse=zeros(0,5);
 for v=10:3:44.5
  for xr=[116,121,126],q=[v,xr,0];score=objective(q,direction);coarse(end+1,:)=[score,q,direction];end
 end
 [~,ix]=sort(coarse(:,1));candidates=[candidates;coarse(ix(1:3),:)]; %#ok<AGROW>
 if matched,candidates(end+1,:)=[objective(railseed,direction),railseed,direction];end
end
best=Inf;
for k=1:size(candidates,1)
 q=candidates(k,2:4);direction=candidates(k,5);fun=@(u)objective(u,direction);
 [qnew,score,iterations]=kw51_nelder_mead(fun,q,[8,114,-.8],[44.44,126,.8],220);
 if score<objective(q,direction),q=qnew;end
 [score,loads,kin]=objective(q,direction);
 if score<best
  best=score;F=struct('kinematics',kin,'rail_position_m',q(2),'loads_N',loads*1e5, ...
   'objective',score,'rail_geometry_match',matched,'optimizer_iterations',iterations);
 end
end
 function [score,loads,kin]=objective(q,direction)
  travel=direction*((direction<0)*116.316-q(2));disc=q(1)^2+2*q(3)*travel;
  loads=[];kin=[];
  if disc<=0,score=1e9;return;end
  dt=2*travel/(q(1)+sqrt(disc));kin=[q(1)+q(3)*dt,E.rail_first_time_s+dt,direction,q(3)];
  H=kw51_basis(B,E,t,kin,refs);[loads,score]=kw51_load_fit(H,y,W);
  if matched,score=score+mean(((q(1)*tau+.5*q(3)*tau.^2-desiro)/.35).^2);end
  if min(kin(1)+kin(4)*(t-kin(2)))<5,score=score+100;end
  if kin(1)>44.44,score=score+100*(kin(1)-44.44)^2;end
 end
end
