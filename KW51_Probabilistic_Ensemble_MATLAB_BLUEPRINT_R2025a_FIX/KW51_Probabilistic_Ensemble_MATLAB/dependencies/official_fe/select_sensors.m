function select_sensors

% Definition of the output selection matrices for the monitoring on railway 
% bridge KW51 in Leuven, Belgium. 
% 
% The sensor layout corresponds to the monitoring described in the 
% following reference:
%
% [1] K. Maes and G. Lombaert. Monitoring railway bridge KW51 before, 
%     during, and after retrofitting. ASCE Journal of Bridge Engineering, 
%     26(3):1-19, 2021. https://doi.org/10.1061/(ASCE)BE.1943-5592.0001668.
%
% The monitoring data obtained from the first 15 months of the monitoring
% are available on Zenodo: https://doi.org/10.5281/zenodo.3745914).
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
%
% # NOTES
% The results are stored in file named 'sensor_definition.mat', which 
% contains the following variables:
% Ss        Selection matrix for the strain measurements.
%           (strain = Ss*U, with U the displacement vector).
%           <12 x ndof> double, with ndof the number of degrees of freedom 
%           in the model.
% Sd        Selection matrix for the displacement measurements. 
%           (displacement = Sd*U).
%           <2 x ndof> double
% Sa        Selection matrix for the acceleration measurements.
%           (displacement = Sa*A, with A the acceleration vector).
%           <12 x ndof> double
% labels_s  Labels of the strain measurements
%           <12 x 1> cell
% labels_d  Labels of the displacement measurements
%           <2 x 1> cell
% labels_a  Labels of the acceleration measurements
%           <12 x 1> cell

%% DEFINITION OF SENSOR LABELS ACCORDING TO [1]

labels_s={'sgBD1011A';'sgBD1415A';'sgBD1718A';'sgBD1718C';'sgBD2324AB';...
    'sgBD2324AT';'sgBD2324C';'sgBD2728A';'sgDI20ALB';'sgDI20ALL';...
    'sgDI23ALB';'sgDI23ALL'};
labels_d={'dABLA';'dABLC'};
labels_a={'aBD11Az','aBD17Ay','aBD17Az','aBD17Cz','aBD23Ay','aBD23Az',...
    'aAR0910Ay','aAR0910Cy','aAR1516Ay','aAR1516Cy','aAR2122Ay','aAR2122Cy'};

%% LOAD FINITE ELEMENT MODEL INFORMATION

% ## Model run with default input parameters (does not affect the results)
[~,~,~,Nodes,Types,Sections,Materials,Elements,DOF]=...
    model_kw51(1,1,0,1,1,0,1,0,[]);

% ## Reference numbers nodes on diagonals A at 0.8 m above the  
% bottom end of the rectangular tube 
% (111XXX - Leaning towards Leuven)
% (112XXX - Leaning towards Brussels)
ndlb=Nodes(floor(Nodes(:,1)/1000)==111|floor(Nodes(:,1)/1000)==112,1);    
% ## Reference numbers nodes on diagonals A at 0.8 m below the  
% top end of the rectangular tube 
% (113XXX - Leaning towards Leuven)
% (114XXX - Leaning towards Brussels)
ndub=Nodes(floor(Nodes(:,1)/1000)==113|floor(Nodes(:,1)/1000)==114,1);    

%% DEFINITION OF STRAIN MEASUREMENTS

% Initialize selection matrix
Ss=zeros(12,length(DOF));

% Material properties
E=Materials(1,2); % Young's modulus of steel

% Section properties main girders
Ag=Sections(1,2);       % Cross section main girder
Iyg=Sections(1,6);      % Bending moment of intertia around local y axis
Izg=Sections(1,7);      % Bending moment of intertia around local z axis
yg=-0.1919;             % Distance between COG and sensor along y in LCS, 
                        % for sensors sgBD1011A, sgBD1415A, sgBD1718A, 
                        % sgBD1718C, sgBD2324C and sgBD2728A
zg=-0.150;              % Distance between COG and sensor along z in LCS, 
                        % for sensors sgBD1011A, sgBD1415A, sgBD1718A, 
                        % sgBD1718C, sgBD2324C and sgBD2728A
ygB=-0.0969;            % Distance between COG and sensor along y in LCS, 
                        % for sensors sgBD2324AB 
zgB=-0.0125;            % Distance between COG and sensor along z in LCS, 
                        % for sensors sgBD2324AB 
ygT=0.3731;             % Distance between COG and sensor along y in LCS, 
                        % for sensors sgBD2324AT 
zgT=-0.0125;            % Distance between COG and sensor along z in LCS, 
                        % for sensors sgBD2324AT 
% Section properties diagonals
Ad=Sections(3,2);       % Cross section diagonal member
Izd=Sections(3,7);      % Bending moment of intertia around local z axis
yd=0.1685;              % Distance between COG and sensor along y in LCS

% Select Nodes on beam elements (exclude plate elements)
eSel=(Elements(:,2)~=3); 
Elements_beam=Elements(eSel,1:end-1);
% List all nodes associated to beam elements
Nodes_beam=unique(reshape(Elements_beam(:,5:6),numel(Elements_beam(:,5:6)),1));
% List all corresponding DOFs
DOF_beam=[Nodes_beam+0.01;Nodes_beam+0.02;Nodes_beam+0.03;
          Nodes_beam+0.04;Nodes_beam+0.05;Nodes_beam+0.06];
% Compute the member forces for unit nodal displacements
[L,I]=selectdof(DOF,DOF_beam); L=L.';
[ForcesLCS,~]=elemforces(Nodes,Elements_beam,Types,Sections,...
    Materials,DOF,L);

% 1 - sgBD1011A
% Select nodes of the corresponding element
iSel1=(floor(Nodes(:,1)/1000)==101)&... % Main girder A
    (round(Nodes(:,2)*10000)==round((34.558-3.933/2)*10000));
iSel2=(floor(Nodes(:,1)/1000)==101)&... % Main girder A
    (round(Nodes(:,2)*10000)==round((34.558)*10000));
Node1=Nodes(iSel1,1); Node2=Nodes(iSel2,1);
% Plot nodes 1 and 2 
plot3(Nodes(iSel1==1,2),Nodes(iSel1==1,3),Nodes(iSel1==1,4),'ko',...
    'MarkerFaceColor','c')
plot3(Nodes(iSel2==1,2),Nodes(iSel2==1,3),Nodes(iSel2==1,4),'ko',...
    'MarkerFaceColor','m')
% Find the corresponding element
eSel=find((Elements_beam(:,5)==Node1)&(Elements_beam(:,6)==Node2)); 
% Find the member forces due to unit nodal displacements
ForcesLCS_=squeeze(ForcesLCS(eSel,:,:));
% Select member forces that contribute to the measured strain (in Node 1)
N=ForcesLCS_(1,:);   % Axial force, NODE 1
Myy=ForcesLCS_(5,:); % Bending moment LCS yy, NODE 1
Mzz=ForcesLCS_(6,:); % Bending moment LCS zz, NODE 1
Ss(1,I)=N/(E*Ag)+Mzz*yg/(E*Izg)+Myy*zg/(E*Iyg);
clear iSel1 iSel2 Node1 Node2 eSel ForcesLCS_ N Myy Mzz

% 2 - sgBD1415A
% Select nodes of the corresponding element
iSel1=(floor(Nodes(:,1)/1000)==101)&...  % Main girder A
    (round(Nodes(:,2)*10000)==round((46.358+3.933/2)*10000));
iSel2=(floor(Nodes(:,1)/1000)==101)&... % Main girder A
    (round(Nodes(:,2)*10000)==round((46.358+3.933)*10000));
Node1=Nodes(iSel1,1); Node2=Nodes(iSel2,1);
% Plot nodes 1 and 2 
plot3(Nodes(iSel1==1,2),Nodes(iSel1==1,3),Nodes(iSel1==1,4),'ko',...
    'MarkerFaceColor','c')
plot3(Nodes(iSel2==1,2),Nodes(iSel2==1,3),Nodes(iSel2==1,4),'ko',...
    'MarkerFaceColor','m')
% Find the corresponding element
eSel=find((Elements_beam(:,5)==Node1)&(Elements_beam(:,6)==Node2)); 
% Find the member forces due to unit nodal displacements
ForcesLCS_=squeeze(ForcesLCS(eSel,:,:));
% Select member forces that contribute to the measured strain (in Node 1)
N=ForcesLCS_(1,:);   % Axial force, NODE 1
Myy=ForcesLCS_(5,:); % Bending moment LCS yy, NODE 1
Mzz=ForcesLCS_(6,:); % Bending moment LCS zz, NODE 1
Ss(2,I)=N/(E*Ag)+Mzz*yg/(E*Izg)+Myy*zg/(E*Iyg);
clear iSel1 iSel2 Node1 Node2 eSel ForcesLCS_ N Myy Mzz

% 3 - sgBD1718A
% Select nodes of the corresponding element
iSel1=(floor(Nodes(:,1)/1000)==101)&...  % Main girder A
    (round(Nodes(:,2)*10000)==round((58.158+3.933/2)*10000));
iSel2=(floor(Nodes(:,1)/1000)==101)&... % Main girder A
    (round(Nodes(:,2)*10000)==round((58.158+3.933/2+0.3)*10000));
Node1=Nodes(iSel1,1); Node2=Nodes(iSel2,1);
% Plot nodes 1 and 2 
plot3(Nodes(iSel1==1,2),Nodes(iSel1==1,3),Nodes(iSel1==1,4),'ko',...
    'MarkerFaceColor','c')
plot3(Nodes(iSel2==1,2),Nodes(iSel2==1,3),Nodes(iSel2==1,4),'ko',...
    'MarkerFaceColor','m')
% Find the corresponding element
eSel=find((Elements_beam(:,5)==Node1)&(Elements_beam(:,6)==Node2)); 
% Find the member forces due to unit nodal displacements
ForcesLCS_=squeeze(ForcesLCS(eSel,:,:));
% Select member forces that contribute to the measured strain (in Node 1)
N=ForcesLCS_(1,:);   % Axial force, NODE 1
Myy=ForcesLCS_(5,:); % Bending moment LCS yy, NODE 1
Mzz=ForcesLCS_(6,:); % Bending moment LCS zz, NODE 1
Ss(3,I)=N/(E*Ag)+Mzz*yg/(E*Izg)+Myy*zg/(E*Iyg);
clear iSel1 iSel2 Node1 Node2 eSel ForcesLCS_ N Myy Mzz

% 4 - sgBD1718C
% Select nodes of the corresponding element
iSel1=(floor(Nodes(:,1)/1000)==102)&... % Main girder C
    (round(Nodes(:,2)*10000)==round((58.158+3.933/2+0.3)*10000));
iSel2=(floor(Nodes(:,1)/1000)==102)&... % Main girder C
    (round(Nodes(:,2)*10000)==round((58.158+3.933)*10000));
Node1=Nodes(iSel1,1); Node2=Nodes(iSel2,1);
% Plot nodes 1 and 2 
plot3(Nodes(iSel1==1,2),Nodes(iSel1==1,3),Nodes(iSel1==1,4),'ko',...
    'MarkerFaceColor','c')
plot3(Nodes(iSel2==1,2),Nodes(iSel2==1,3),Nodes(iSel2==1,4),'ko',...
    'MarkerFaceColor','m')
% Find the corresponding element
eSel=find((Elements_beam(:,5)==Node1)&(Elements_beam(:,6)==Node2)); 
% Find the member forces due to unit nodal displacements
ForcesLCS_=squeeze(ForcesLCS(eSel,:,:));
% Select member forces that contribute to the measured strain (in Node 1)
N=ForcesLCS_(1,:);   % Axial force, NODE 1
Myy=ForcesLCS_(5,:); % Bending moment LCS yy, NODE 1
Mzz=ForcesLCS_(6,:); % LCS zz, NODE 1
Ss(4,I)=N/(E*Ag)+Mzz*yg/(E*Izg)-Myy*zg/(E*Iyg);
clear iSel1 iSel2 Node1 Node2 eSel ForcesLCS_ N Myy Mzz

% 5 - sgBD2324AB
% Select nodes of the corresponding element
iSel1=(floor(Nodes(:,1)/1000)==101)&...  % Main girder A
    (round(Nodes(:,2)*10000)==round((81.758+3.933/2)*10000));
iSel2=(floor(Nodes(:,1)/1000)==101)&... % Main girder A
    (round(Nodes(:,2)*10000)==round((81.758+3.933)*10000));
Node1=Nodes(iSel1,1); Node2=Nodes(iSel2,1);
% Plot nodes 1 and 2 
plot3(Nodes(iSel1==1,2),Nodes(iSel1==1,3),Nodes(iSel1==1,4),'ko',...
    'MarkerFaceColor','c')
plot3(Nodes(iSel2==1,2),Nodes(iSel2==1,3),Nodes(iSel2==1,4),'ko',...
    'MarkerFaceColor','m')
% Find the corresponding element
eSel=find((Elements_beam(:,5)==Node1)&(Elements_beam(:,6)==Node2)); 
% Find the member forces due to unit nodal displacements
ForcesLCS_=squeeze(ForcesLCS(eSel,:,:));
% Select member forces that contribute to the measured strain (in Node 1)
N=ForcesLCS_(1,:);   % Axial force, NODE 1
Myy=ForcesLCS_(5,:); % Bending moment LCS yy, NODE 1
Mzz=ForcesLCS_(6,:); % Bending moment LCS zz, NODE 1
Ss(5,I)=N/(E*Ag)+Mzz*ygB/(E*Izg)+Myy*zgB/(E*Iyg);

% 6 - sgBD2324AT
Ss(6,I)=N/(E*Ag)+Mzz*ygT/(E*Izg)+Myy*zgT/(E*Iyg);
clear iSel1 iSel2 Node1 Node2 eSel ForcesLCS_ N Myy Mzz

% 7 - sgBD2324C
iSel1=(floor(Nodes(:,1)/1000)==102)&...  % Main girder C
    (round(Nodes(:,2)*10000)==round((81.758+3.933/2)*10000));
iSel2=(floor(Nodes(:,1)/1000)==102)&... % Main girder C
    (round(Nodes(:,2)*10000)==round((81.758+3.933)*10000));
Node1=Nodes(iSel1,1); Node2=Nodes(iSel2,1);
% Plot nodes 1 and 2 
plot3(Nodes(iSel1==1,2),Nodes(iSel1==1,3),Nodes(iSel1==1,4),'ko',...
    'MarkerFaceColor','c')
plot3(Nodes(iSel2==1,2),Nodes(iSel2==1,3),Nodes(iSel2==1,4),'ko',...
    'MarkerFaceColor','m')
% Find the corresponding element
eSel=find((Elements_beam(:,5)==Node1)&(Elements_beam(:,6)==Node2)); 
% Find the member forces due to unit nodal displacements
ForcesLCS_=squeeze(ForcesLCS(eSel,:,:));
% Select member forces that contribute to the measured strain (in Node 1)
N=ForcesLCS_(1,:);   % Axial force, NODE 1
Myy=ForcesLCS_(5,:); % Bending moment LCS yy, NODE 1
Mzz=ForcesLCS_(6,:); % Bending moment LCS zz, NODE 1
Ss(7,I)=N/(E*Ag)+Mzz*yg/(E*Izg)-Myy*zg/(E*Iyg);
clear iSel1 iSel2 Node1 Node2 eSel ForcesLCS_ N Myy Mzz

% 8 - sgBD2728A
iSel1=(floor(Nodes(:,1)/1000)==101)&... % Main girder A
    (round(Nodes(:,2)*10000)==round((97.391+3.834/2)*10000));
iSel2=(floor(Nodes(:,1)/1000)==101)&... % Main girder A
    (round(Nodes(:,2)*10000)==round((97.391+3.834)*10000));
Node1=Nodes(iSel1,1); Node2=Nodes(iSel2,1);
% Plot nodes 1 and 2 
plot3(Nodes(iSel1==1,2),Nodes(iSel1==1,3),Nodes(iSel1==1,4),'ko',...
    'MarkerFaceColor','c')
plot3(Nodes(iSel2==1,2),Nodes(iSel2==1,3),Nodes(iSel2==1,4),'ko',...
    'MarkerFaceColor','m')
% Find the corresponding element
eSel=find((Elements_beam(:,5)==Node1)&(Elements_beam(:,6)==Node2)); 
% Find the member forces due to unit nodal displacements
ForcesLCS_=squeeze(ForcesLCS(eSel,:,:));
% Select member forces that contribute to the measured strain (in Node 1)
N=ForcesLCS_(1,:);   % Axial force, NODE 1
Myy=ForcesLCS_(5,:); % Bending moment LCS yy, NODE 1
Mzz=ForcesLCS_(6,:); % Bending moment LCS zz, NODE 1
Ss(8,I)=N/(E*Ag)+Mzz*yg/(E*Izg)+Myy*zg/(E*Iyg);
clear iSel1 iSel2 Node1 Node2 eSel ForcesLCS_ N Myy Mzz

% 9 - sgDI20ALB
Node1=ndlb(6); Node2=ndub(6);
% Plot nodes 1 and 2 
plot3(Nodes(Nodes(:,1)==Node1,2),Nodes(Nodes(:,1)==Node1,3),...
    Nodes(Nodes(:,1)==Node1,4),'ko','MarkerFaceColor','c')
plot3(Nodes(Nodes(:,1)==Node2,2),Nodes(Nodes(:,1)==Node2,3),...
    Nodes(Nodes(:,1)==Node2,4),'ko','MarkerFaceColor','m')
% Find the corresponding element
eSel=find((Elements_beam(:,5)==Node1)&(Elements_beam(:,6)==Node2)); 
% Find the member forces due to unit nodal displacements
ForcesLCS_=squeeze(ForcesLCS(eSel,:,:));
N=ForcesLCS_(1,:);
Mzz=ForcesLCS_(6,:);
Ss(9,I)=N/(E*Ad)+Mzz*yd/(E*Izd);

% 10 - sgDI20ALL
Ss(10,I)=N/(E*Ad)-Mzz*yd/(E*Izd);
clear Node1 Node2 eSel ForcesLCS_ N Mzz

% 11 - sgDI23ALB
Node1=ndlb(7); Node2=ndub(7);
% Plot nodes 1 and 2 
plot3(Nodes(Nodes(:,1)==Node1,2),Nodes(Nodes(:,1)==Node1,3),...
    Nodes(Nodes(:,1)==Node1,4),'ko','MarkerFaceColor','c')
plot3(Nodes(Nodes(:,1)==Node2,2),Nodes(Nodes(:,1)==Node2,3),...
    Nodes(Nodes(:,1)==Node2,4),'ko','MarkerFaceColor','m')
% Find the corresponding element
eSel=find((Elements_beam(:,5)==Node1)&(Elements_beam(:,6)==Node2)); 
% Find the member forces due to unit nodal displacements
ForcesLCS_=squeeze(ForcesLCS(eSel,:,:));
N=ForcesLCS_(1,:);   % Axial force, NODE 1
Mzz=ForcesLCS_(6,:); % Bending moment LCS zz, NODE 1
Ss(11,I)=N/(E*Ad)+Mzz*yd/(E*Izd);

% 12 - sgDI23ALL
Ss(12,I)=N/(E*Ad)-Mzz*yd/(E*Izd);

%% DEFINITION OF DISPLACEMENT MEASUREMENTS

Sd=selectdof(DOF,[3.01;4.01]);

%% DEFINITION OF ACCELERATION MEASUREMENTS

% Initialize selection matrix
Sa=zeros(12,length(DOF));

% 1 - aBD11Az
iSel=(floor(Nodes(:,1)/1000)==101)&... % Main girder A
    (round(Nodes(:,2)*10000)==round(34.558*10000));
hold on
plot3(Nodes(iSel,2),Nodes(iSel,3),Nodes(iSel,4),'ko',...
    'MarkerFaceColor','r')
Sa(1,:)=selectdof(DOF,Nodes(iSel,1)+0.03);
clear iSel

% 2 - aBD17Ay
iSel=(floor(Nodes(:,1)/1000)==101)&... % Main girder A
    (round(Nodes(:,2)*10000)==round(58.158*10000));
plot3(Nodes(iSel,2),Nodes(iSel,3),Nodes(iSel,4),'ko',...
    'MarkerFaceColor','r')
Sa(2,:)=selectdof(DOF,Nodes(iSel,1)+0.02)...
    +selectdof(DOF,Nodes(iSel,1)+0.04)*(271.94-150)/1000;

% 3 - aBD17Az
Sa(3,:)=selectdof(DOF,Nodes(iSel,1)+0.03);
clear iSel

% 4 - aBD17Cz
iSel=(floor(Nodes(:,1)/1000)==102)&... % Main girder C
    (round(Nodes(:,2)*10000)==round(58.158*10000));
plot3(Nodes(iSel,2),Nodes(iSel,3),Nodes(iSel,4),'ko',...
    'MarkerFaceColor','r')
Sa(4,:)=selectdof(DOF,Nodes(iSel,1)+0.03);
clear iSel

% 5 - aBD23Ay
iSel=(floor(Nodes(:,1)/1000)==101)&... % Main girder A
    (round(Nodes(:,2)*10000)==round(81.758*10000));
plot3(Nodes(iSel,2),Nodes(iSel,3),Nodes(iSel,4),'ko',...
    'MarkerFaceColor','r')
Sa(5,:)=selectdof(DOF,Nodes(iSel,1)+0.02)...
    +selectdof(DOF,Nodes(iSel,1)+0.04)*(271.94-150)/1000;

% 6 - aBD23Az
Sa(6,:)=selectdof(DOF,Nodes(iSel,1)+0.03);
clear iSel

% 7 - aAR0910Ay
iSel=(floor(Nodes(:,1)/1000)==103)&... % Arch A
    (round(Nodes(:,2)*10000)==round(28.658*10000));
plot3(Nodes(iSel,2),Nodes(iSel,3),Nodes(iSel,4),'ko',...
    'MarkerFaceColor','r')
Sa(7,:)=selectdof(DOF,Nodes(iSel,1)+0.02)...
    +selectdof(DOF,Nodes(iSel,1)+0.04)*(271.94-150)/1000;
clear iSel

% 8 - aAR0910Cy
iSel=(floor(Nodes(:,1)/1000)==104)&... % Arch C
    (round(Nodes(:,2)*10000)==round(28.658*10000));
plot3(Nodes(iSel,2),Nodes(iSel,3),Nodes(iSel,4),'ko',...
    'MarkerFaceColor','r')
Sa(8,:)=selectdof(DOF,Nodes(iSel,1)+0.02)...
    +selectdof(DOF,Nodes(iSel,1)+0.04)*(271.94-150)/1000;
clear iSel

% 9 - aAR1516Ay
iSel=(floor(Nodes(:,1)/1000)==103)&... % Arch A
    (round(Nodes(:,2)*10000)==round(52.258*10000));
plot3(Nodes(iSel,2),Nodes(iSel,3),Nodes(iSel,4),'ko',...
    'MarkerFaceColor','r')
Sa(9,:)=selectdof(DOF,Nodes(iSel,1)+0.02)...
    +selectdof(DOF,Nodes(iSel,1)+0.04)*(271.94-150)/1000;
clear iSel

% 10 - aAR1516Cy
iSel=(floor(Nodes(:,1)/1000)==104)&... % Arch C
    (round(Nodes(:,2)*10000)==round(52.258*10000));
plot3(Nodes(iSel,2),Nodes(iSel,3),Nodes(iSel,4),'ko',...
    'MarkerFaceColor','r')
Sa(10,:)=selectdof(DOF,Nodes(iSel,1)+0.02)...
    +selectdof(DOF,Nodes(iSel,1)+0.04)*(271.94-150)/1000;
clear iSel

% 11 - aAR2122Ay
iSel=(floor(Nodes(:,1)/1000)==103)&... % Arch A
    (round(Nodes(:,2)*10000)==round(75.858*10000));
plot3(Nodes(iSel,2),Nodes(iSel,3),Nodes(iSel,4),'ko',...
    'MarkerFaceColor','r')
Sa(11,:)=selectdof(DOF,Nodes(iSel,1)+0.02)...
    +selectdof(DOF,Nodes(iSel,1)+0.04)*(271.94-150)/1000;
clear iSel

% 12 - aAR2122Cy
iSel=(floor(Nodes(:,1)/1000)==104)&... % Arch C
    (round(Nodes(:,2)*10000)==round(75.858*10000));
plot3(Nodes(iSel,2),Nodes(iSel,3),Nodes(iSel,4),'ko',...
    'MarkerFaceColor','r')
Sa(12,:)=selectdof(DOF,Nodes(iSel,1)+0.02)...
    +selectdof(DOF,Nodes(iSel,1)+0.04)*(271.94-150)/1000;
clear iSel

%% Save the results
save('sensor_definition.mat',...
    'Ss','Sd','Sa','labels_s','labels_d','labels_a')

end