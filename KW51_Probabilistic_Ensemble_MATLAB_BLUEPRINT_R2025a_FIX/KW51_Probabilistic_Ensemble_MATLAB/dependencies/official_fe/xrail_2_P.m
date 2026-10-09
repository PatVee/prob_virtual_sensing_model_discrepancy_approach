function P=xrail_2_P(track,xrail,AL,va)

% Construct the equivalent load vectors for given axle loads at given 
% locations and for given train speeds, taking into account both
% gravitational and centripetal forces. Multiple time steps can be
% considered. The forces corresponding to the different axles for a given
% time step are combined in a single load vector (column of P).
%
% The load distribution is described in the following reference:
%
% [1] K. Maes and G. Lombaert. Validation of virtual sensing for the recon-
%     struction of stresses in a railway bridge using field data of the  
%     KW51 bridge. Mechanical Systems and Signal Processing, 190(110142):
%     1-16, 2023. https://doi.org/10.1016/j.ymssp.2023.110142.
%
% October 2025 
%
% # AUTHORS
% Kristof Maes, Menno van de Velde and Geert Lombaert
% KU Leuven, Department of Civil Engineering, Structural Mechanics Section
%
% # UNITS 
% [kg, m, s]
%
% # INPUTS
% track     Track number, 'A' or 'B' (see overview.pdf)
%           <1 x 1> char
% xrail     Axle locations (curvilinear measurement) with respect to the
%           onset of the bridge at x=0.
%           <N x Na> double, with N the number of time steps and Na the
%           number of axles.
% AL        Axle load vector (axle loads in Newton)
%           <1 x Na> double
% va        Train speed 
%           <N x 1> double
%
% # OUTPUTS
% P         Load vectors (columns of P) for all time steps (rows of P)
%           <Ndof x N> double

%% DEFAULT INPUT ARGUMENTS

if nargin==0
    track='B';
    xrail=(-2:118).';
    AL=1;
    va=ones(size(xrail));
    % Illustrate force distribution
    force_distr;
end

%% ROBUSTNESS CHECKS

% Track: must be 'A' or 'B' 
if ~ismember(track, {'A','B'})
    error('Invalid track ''%s''. Must be ''A'' or ''B''.', track);
end

% xrail
if ~isnumeric(xrail) || ~ismatrix(xrail)
    error('xrail must be a real numeric 2-D array of size N×Na.');
end
if any(~isfinite(xrail(:)))
    error('xrail contains NaN or Inf.');
end
[N, Na] = size(xrail);

% AL
if ~isnumeric(AL) || ~isreal(AL)
    error('AL must be a real numeric vector (forces in N).');
end
if numel(AL) ~= Na
    error('AL must have Na elements (got %d, expected %d).', numel(AL), Na);
end

% va
if ~isnumeric(va) || ~isreal(va)
    error('va must be a real numeric vector (m/s).');
end
if numel(va) ~= N
    error('va must have N elements (got %d, expected %d).', numel(va), N);
end

% ## Convert curvilinear axle location to xP and yP in the global coordi-
% nate system (GCS) of the model (see overview.pdf).
[xP,yP]=xrail_2_xPyP(track,xrail);

% ## Compose the load vector, in accordance to the load distribution 
% described in [1]. 

P=[]; % Initialize
% Curve radius of the tracks
if strcmp(track,'B'), Rad=1121;
elseif strcmp(track,'A'),  Rad=1125;
end
% Distance between the center of the ballast to the mass COG of the train.
H=2.5;  
% Compose the load vector
for indx=1:size(xrail,1)    % Loop over all time steps
    for inda=1:length(AL)   % Loop over all axles
        % Vertical load along z (gravitational force)
        Pz=-AL(inda)*force_distr(xP(indx,inda),yP(indx,inda),0.03,0);
        % Lateral load along y (centripetal force)
        Py=+AL(inda)/9.81*va(indx)^2/Rad*...
            force_distr(xP(indx,inda),yP(indx,inda),0.02,0);
        % Moment around x (centripetal force)
        Mx=-AL(inda)/9.81*va(indx)^2/Rad*H*...
            force_distr(xP(indx,inda),yP(indx,inda),0.06,0);
        % Combine the three contributions
        if inda==1
            Pc=Pz+Py+Mx;
        else
            Pc=Pc+Pz+Py+Mx;
        end
    end
    P=[P,Pc];
    clear Pc
end

end

%% AUXILIARY FUNCTIONS

% ## ====================================================================##
function [xP,yP]=xrail_2_xPyP(track,xrail)

% Convert curvilinear axle location to xP and yP in the global coordinate
% system (GCS) of the model (see overview.pdf).
%
% October 2025 
%
% # AUTHORS
% Kristof Maes, Menno van de Velde and Geert Lombaert
% KU Leuven, Department of Civil Engineering, Structural Mechanics Section
%
% # UNITS 
% [kg, m, s]
%
% # INPUTS
% track     Track number, 'A' or 'B'
%           <1 x 1> char
% xrail     Axle locations (curvilinear measurement) with respect to the
%           onset of the bridge at x=0.
%           <N x Na> double, with N the number of time steps and Na the
%           number of axles.
%
% # OUTPUTS
% xP        Axle location along x in GCS
%           <N x Na> double
% yP        Axle location along y in GCS
%           <N x Na> double
%
% # NOTES
%
% - This function can be applied for multiple (Na) axle locations under the 
%   condition that they are on the same track.

% Curve fit center line of track B, based on blueprints
xBr=[0,58.158,115.658];
yBr=3.24+[0,1.56,0];
PB=polyfit(xBr,yBr,2);
% Curve fit center line of track A, based on blueprints
xAr=xBr; yAr=7.24+[0,1.57,0];
PA=polyfit(xAr,yAr,2);

xP=xrail; % Curvature of the rail is disregarded
if strcmp(track,'B')
    yP=polyval(PB,xP);
elseif strcmp(track,'A')
    yP=polyval(PA,xP);
end

end

% ## ====================================================================##
function P=force_distr(xP,yP,excdir,illustr)

% Construct the equivalent load vector for a unit load applied at location
% (xP,yP) in the global coordinate system (GCS) of the model (see 
% overview.pdf).
%
% October 2025 
%
% # AUTHORS
% Kristof Maes, Menno van de Velde and Geert Lombaert
% KU Leuven, Department of Civil Engineering, Structural Mechanics Section
%
% # UNITS 
% [kg, m, s]
%
% # INPUTS
% xP        Axle location along x in GCS
%           <1 x 1> double
% yP        Axle location along y in GCS
%           <1 x 1> double
% excdir    Excitation direction in Stabil convention
%           0.01 = point load along x in GCS
%           0.02 = point load along y in GCS
%           0.03 = point load along z in GCS
%           0.04 = bending moment around x in GCS
%           0.05 = bending moment around y in GCS
%           0.06 = bending moment around z in GCS
%           <1 x 1> double
% illustr   Load distribution plot indicator.
%           <1 x 1> double, 0 or 1
%           1 - plot illustration of the load distribution.
%           0 - do not plot illustration.
%
% # OUTPUTS
% P         Load vector
%           <Ndof x 1> double, with Ndof the number of degrees of freedom
%           in the model.
%
% # NOTES
%
% - This function can only be applied for a single force location (xP,yP).
% 
% - To obtain the illustration in force_distribution.pdf, run this function
%   without input arguments.

% Default input arguments
if nargin==0
    xP=55; yP=5.8;
    excdir=0.03;
    illustr=1;
end

% Influence zone of a load along the x-direction (longitudinal).
% Assumption: Load is taken by 7 sleepers at an intermediate distance of 
% 0.6 m: - - - + - - - => 0.6*6/2 m 
% + distribution in the ballast under an angle of 45 degrees => + 0.6 m
infl_x=0.6*4; 

% Influence zone of a load along the y-direction (lateral).
% Distance between two rails = 1.435 m => 1.435/2 m
% + distribution at 45 degrees => + 0.6 m
infl_y=1.435/2+0.6;

% Compute or load the areas assigned to each node of the ballast layer.
if ~isfile('plate_areas.mat')
    plate_areas; % Compute plate areas and store them in plate_areas.mat
end

% Load plate_areas.mat
load('plate_areas.mat','NodesSel','yL','yR','xT','xB','Nodes','DOF')

if illustr==1

    % Curve fit center line of track B, based on blueprints
    xBr=[0,58.158,115.658];
    yBr=3.24+[0,1.56,0];
    PB=polyfit(xBr,yBr,2);
    % Curve fit center line of track A, based on blueprints
    xAr=xBr; yAr=7.24+[0,1.57,0];
    PA=polyfit(xAr,yAr,2);
    xr=0:0.01:115;
    yrA=polyval(PA,xr);
    yrB=polyval(PB,xr);

    figure % Open new figure
    for ind1=1:length(NodesSel)
        iSel=Nodes(:,1)==NodesSel(ind1);
        plot(Nodes(iSel,2),Nodes(iSel,3),'b.','markersize',5);
        hold on
    end
    if nargin==0 % Adjust y-location of the load (on center rail)
        [~,imin]=min(abs(xr-xP)); yP=yrB(imin); clear imin
    end
    axis equal
    xlabel('x GCS [m]')
    ylabel('y GCS [m]')
    xlim(xP+[-4.6,4.6])
    ylim(yP+[-3.6,3.6])
    yL0=yL; yR0=yR;
    xB0=xB; xT0=xT;
    
end

% Boundaries of the load influence zone
yLI=yP+infl_y; % Left
yRI=yP-infl_y; % Right
xTI=xP+infl_x; % Top 
xBI=xP-infl_x; % Bottom

for ind1=1:length(NodesSel)
    if yL(ind1)>yLI, yL(ind1)=yLI; end  
    if yR(ind1)<yRI, yR(ind1)=yRI; end
    if xT(ind1)>xTI, xT(ind1)=xTI; end
    if xB(ind1)<xBI, xB(ind1)=xBI; end    
end

% Select nodes of the ballast layer where the load is transferred
iSel=(yL>yR)&(xT>xB); 
% Corresponding node numbers
Nodes_p=NodesSel(iSel,:); 
% Corresponding boundaries of the assigned rectangular areas
yL=yL(iSel); yR=yR(iSel); xT=xT(iSel); xB=xB(iSel);
% Surface areas
Ap=(yL-yR).*(xT-xB);
% Scaling factor when load is not fully transferred to the bridge
% (axle enters or leaves the bridge)
if xBI<0
    scalfac=(xTI-0)/(xTI-xBI);
elseif xTI>116.316
    scalfac=(116.316-xBI)/(xTI-xBI);
else 
    scalfac=1;
end
% Load vector
P=nodalvalues(DOF,Nodes_p+excdir,scalfac*Ap/sum(Ap));

if illustr==1

    yL0=yL0(iSel); yR0=yR0(iSel); xB0=xB0(iSel); xT0=xT0(iSel);

    p0=plot(xP+[-infl_x,+infl_x,+infl_x,-infl_x,-infl_x,+infl_x],...
        yP+[-infl_y,-infl_y,+infl_y,+infl_y,-infl_y,-infl_y],'-',...
        'color',[0,1,0],'linewidth',5);
    for ind1=1:length(Nodes_p)
        p1=plot([xB0(ind1),xT0(ind1),xT0(ind1),xB0(ind1),xB0(ind1)],...
            [yR0(ind1),yR0(ind1),yL0(ind1),yL0(ind1),yR0(ind1)],'b-');
        rectangle('Position',...
            [xB(ind1) yR(ind1) xT(ind1)-xB(ind1) yL(ind1)-yR(ind1)],...
            'EdgeColor','r','FaceColor',[0.8,0.8,0.8]);
        iSel=Nodes(:,1)==Nodes_p(ind1);
        p2=plot(Nodes(iSel,2),Nodes(iSel,3),'b.','markersize',5);
    end
    p3=patch(NaN,NaN,[0.8 0.8 0.8],'EdgeColor','r'); % Dummy patch for legend
    p4=plot(xP,yP,'r.','markersize',20);
    plot(xr,yrA,'k--');
    plot(xr,yrA-1.435/2,'k-');
    plot(xr,yrA+1.435/2,'k-');
    plot(xr,yrB,'k--');
    plot(xr,yrB-1.435/2,'k-');
    plot(xr,yrB+1.435/2,'k-');
    legend([p2,p1,p4,p0,p3],...
        {'FEM Nodes','Area assigned to node','Force location',...
        'Influence area of force','Effectively assigned area'},...
        'NumColumns',2)
    title('Top view bridge deck illustrating the force distribution')

end

end

% ## ====================================================================##
function []=plate_areas

% Calculate the area that is assigned to each node of the ballast layer. 
% This is the basis for the load distribution. 
%
% October 2025 
%
% # AUTHORS
% Kristof Maes, Menno van de Velde and Geert Lombaert
% KU Leuven, Department of Civil Engineering, Structural Mechanics Section
%
% # UNITS 
% [kg, m, s]
%
% # INPUTS
% <None>
%
% # OUTPUTS
% <None>

% Load Nodes and degrees of freedom (DOF) from the FE model
[~,~,~,Nodes,~,~,~,~,DOF]=model_kw51(1,1,NaN,1,1,1,0,0,[]);

% Select all nodes of the ballast layer
% Numbering of nodes of ballast layer:
% 282000 - 2830XX - 2835XX - ... - 2940XX - 2945XX - 2950XX
% XX = 1:78
numSely1=sort([2820:10:2950,2835:10:2945]);
iSel=[];
for ind1=1:length(numSely1)
    iSel=[iSel;find(floor(Nodes(:,1)/100)==numSely1(ind1))];
end
Nodes=Nodes(iSel,:);

% Compute assigned area for all internal nodes 
numSely2=[2830:10:2940,2835:10:2945];
numSelx2=2:76;
NodesSel=zeros(length(numSely2)*length(numSelx2),1);
yL=zeros(size(NodesSel)); yR=zeros(size(NodesSel)); 
xT=zeros(size(NodesSel)); xB=zeros(size(NodesSel));
ind=0;
newfig
for ind1=1:length(numSely2)
    for ind2=1:length(numSelx2)
        ind=ind+1;
        % Select node
        NodesSel(ind)=100*numSely2(ind1)+numSelx2(ind2);
        Node_=NodesSel(ind);
        iSel=find(numSely1==floor(Node_/100));
        % First node to the left
        Node_L=100*numSely1(iSel+1)+numSelx2(ind2);
        % First node to the right
        Node_R=100*numSely1(iSel-1)+numSelx2(ind2);
        % First node to the top
        Node_T=Node_+1;
        % First node to the bottom
        Node_B=Node_-1;        
        Node_=Nodes(Nodes(:,1)==Node_,:);
        yL(ind)=Node_(3)+1/2*(Nodes(Nodes(:,1)==Node_L,3)-Node_(3));
        yR(ind)=Node_(3)+1/2*(Nodes(Nodes(:,1)==Node_R,3)-Node_(3));
        xT(ind)=Node_(2)+1/2*(Nodes(Nodes(:,1)==Node_T,2)-Node_(2));
        xB(ind)=Node_(2)+1/2*(Nodes(Nodes(:,1)==Node_B,2)-Node_(2));
        
    end
end

% Save results
save('plate_areas.mat','NodesSel','yL','yR','xT','xB','Nodes','DOF')

end
% ## ====================================================================##