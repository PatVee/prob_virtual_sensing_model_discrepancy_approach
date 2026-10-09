function [Kr,Mr,T,Nodes,Types,Sections,Materials,Elements,DOF]=model_kw51(Es_rel,hb_rel,kxx,Id_rel,Ad_rel,fixed,plot_elem,addmass,diags_box)

% Finite Element model of railway bridge KW51 in Leuven, Belgium.
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
% Es_rel    Scaling factor for the Young's modulus of steel, Es.
%           <1 x 1> double
%           (default: 1, corresponding to Es = 200 GPa)
% hb_rel    Scaling factor for the thickness of the ballast layer, hb.
%           <1 x 1> double
%           (default: 1, corresponding to hb = 0.6 m)
% kxx       Stiffness of longitudinal springs at the bearings side Leuven.
%           <1 x 1> double
%           (default: 1e13 N/m)
% Id_rel    Scaling factor for the out-of-plane bending moment of inertia  
%           I of the diagonal members at the connections with the arches 
%           and the bridge deck.
%           <1 x 1> double
%           (default: 1, corresponding to the stiffness of the rectangular 
%           tube (I = 3.840e-4 m^4))
% Ad_rel    Scaling factor for the cross section A of the diagonal  members
%           at the connections with the arches and the bridge deck. 
%           <1 x 1> double
%           (default: 1, corresponding to the cross section of the
%           rectangular tube (A = 2.089e-2 m^2))
% fixed     Boundary condition indicator.
%           <1 x 1> double, 0 or 1
%           1 - bridge is pinned at the four bearings (default).
%           0 - roller supports at side Leuven and boundary conditions as
%               indicated in blueprint (KW51 - 13) - GTN-D0-L201-03-015.
% plot_elem Geometry (nodes and elements) plot indicator.
%           <1 x 1> double, 0 or 1
%           1 - plot geometry (nodes and elements, default).
%           0 - do not plot geometry.
% addmass   Indicator for the added mass of the retrofit boxes.
%           <1 x 1> double, 0 or 1
%           0 - do not include the mass of the boxes (default, corresponds 
%               to the situation before the retrofit).
%           1 - add the mass of the boxes (333 kg/box, corresponds to the 
%               situation after the retrofit).
% diags_box List of numbers of the diagonals where a box is modeled, both
%           in terms of added mass and stiffness.
%           10XX - track B, bottom, diag XX (numbering, see overview.pdf)
%           20XX - track B, top, diag XX
%           30XX - track A, bottom, diag XX
%           40XX - track A, top, diag XX
%           (default: empty, [])
%           (example: [1001, 1005, 1006])
%
% # OUTPUTS
% Kr        Stiffness matrix for condensed set of degrees of freedom.
%           <ndof_r x ndof_r> sparse double, with ndof_r the number of
%           reduced degrees of freedom.
% Mr        Mass matrix for condensed set of degrees of freedom.
%           <ndof_r x ndof_r> sparse double, with ndof_r the number of
%           reduced degrees of freedom.
% T         Transformation matrix according to the condensation. 
%           <ndof x ndof_r> sparse double, with ndof the total number of 
%           degrees of freedom and ndof_r the number of reduced degrees of 
%           freedom.
%           The following transformations apply (see also help tconstr.m):
%           Stiffness matrix:    Kr=T.'*K*T 
%           Mass matrix:         Mr=T.'*M*T 
%           Load vector:         Fr=T.'*(F-K*Q0)  (Q0 is zero in this case)
%           Displacement vector: U=T*Ur+Q0
% Nodes     Node definitions          [NodID x y z] 
%           (see Stabil manual)
%           <nnode x 4> double
% Types     Element type definitions  {TypID EltName Option1 ...}
%           (see Stabil manual)
%           <ntype x 2> cell
% Sections  Section definitions       [SecID SecProp1 SecProp2 ...]
%           (see Stabil manual)
%           <nsect x 7> double
% Materials Material definitions      [MatID MatProp1 MatProp2 ...]
%           (see Stabil manual)
%           <nmat x 4> double
% Elements  Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%           (see Stabil manual)
%           <nelem x 7> double
% DOF       Vector with all degrees of freedom in accordance to the Stabil
%           convention (see Stabil manual).
%           <ndof x 1> double
%
% # IMPORTANT NOTES
%
% - Before using the model, please consult 'overview.pdf' for an overview 
%   of the geometry, the definition of the reference coordinate system, the
%   labeling of the bridge sides (South = A, North = C), the numbering of
%   the diagonals (1 - 16), and the bearings (N1 - N4).
%
% - The blueprints of the bridge underlying the model are provided in the 
%   folder 'blueprints', with a description of the blueprints in
%   'list_of_blueprints.pdf'.
%
% - The model is accompanied by two application examples:
%
%   1) 'simulate_trainpassage.m'
%
%       Quasi-static simulation of a train passage, in accordance to the
%       following reference:
%
%       [1] K. Maes and G. Lombaert. Validation of virtual sensing for the 
%           reconstruction of stresses in a railway bridge using field data  
%           of the KW51 bridge. Mechanical Systems and Signal Processing, 
%           190(110142):1-16, 2023. 
%           https://doi.org/10.1016/j.ymssp.2023.110142.
%
%   2) 'modal_analysis.m'
%
%       Computation of the modal characteristics, in accordance to the 
%       following reference:
%
%       [2] K. Maes, L. Van Meerbeeck, E. Reynders, and G. Lombaert.  
%           Validation of vibration-based structural health monitoring on 
%           retrofitted railway bridge KW51. Mechanical Systems and Signal 
%           Processing, 165(108380):1-24, 2022. 
%           https://doi.org/10.1016/j.ymssp.2021.108380.
%
% - The data from the first 15 months of the monitoring on railway bridge 
%   KW51 are available at the following Zenodo repository:
%
%       [3] K. Maes and G. Lombaert. Monitoring data for railway bridge 
%           KW51 in Leuven, Belgium, before, during, and after retrofitting 
%           (Version 1.0), Zenodo, 2020. 
%           https://doi.org/10.5281/zenodo.3745914.
% 
% - A full description of the data available on Zenodo [3] is found in the 
%   following reference:
%  
%       [4] K. Maes and G. Lombaert. Monitoring railway bridge KW51 before, 
%           during, and after retrofitting. ASCE Journal of Bridge Enginee-
%           ring, 26(3):1-19, 2021. 
%           https://doi.org/10.1061/(ASCE)BE.1943-5592.0001668.

%% DEFAULT INPUT ARGUMENTS

if nargin==0, Es_rel=1;     end
if nargin<2,  hb_rel=1;     end
if nargin<3,  kxx=1e13;     end
if nargin<4,  Id_rel=1;     end
if nargin<5,  Ad_rel=1;     end 
if nargin<6,  fixed=1;      end
if nargin<7,  plot_elem=1;  end
if nargin<8,  addmass=0;    end
if nargin<9,  diags_box=[]; end

%% ROBUSTNESS CHECKS

assert(isscalar(Es_rel) && isnumeric(Es_rel) && isfinite(Es_rel) && Es_rel>0, ...
    'Es_rel must be a positive finite scalar.');
assert(isscalar(hb_rel) && isnumeric(hb_rel) && isfinite(hb_rel) && hb_rel>0, ...
    'hb_rel must be a positive finite scalar.');
assert(isscalar(Id_rel) && isnumeric(Id_rel) && isfinite(Id_rel) && Id_rel>0, ...
    'Id_rel must be a positive finite scalar.');
assert(isscalar(Ad_rel) && isnumeric(Ad_rel) && isfinite(Ad_rel) && Ad_rel>0, ...
    'Ad_rel must be a positive finite scalar.');

% Flags ∈ {0,1}
assert(isscalar(fixed)     && ismember(fixed,    [0 1]), 'fixed must be 0 or 1.');
assert(isscalar(plot_elem) && ismember(plot_elem,[0 1]), 'plot_elem must be 0 or 1.');
assert(isscalar(addmass)   && ismember(addmass,  [0 1]), 'addmass must be 0 or 1.');

% kxx
if fixed==0
    assert(isscalar(kxx) && isnumeric(kxx) && isfinite(kxx) && kxx>=0, ...
        'kxx must be a nonnegative finite scalar when fixed==0.');
end

% diags_box: empty or vector of allowed codes
if ~isempty(diags_box)
    assert(isnumeric(diags_box) && isvector(diags_box) && all(isfinite(diags_box)), ...
        'diags_box must be empty or a finite numeric vector.');
    grp = floor(diags_box/100);
    idx = diags_box - grp*100;
    assert(all(ismember(grp,[10 20 30 40])) && all(idx>=1 & idx<=16), ...
        'diags_box entries must be 10xx/20xx/30xx/40xx with 01≤xx≤16.');
end

%% NODES DEFINITION (GEOMETRY ACCORDING TO BLUEPRINTS)

% Nodes=[NodID x y z]
% Reference coordinate system, see overview.pdf

% Ballast thickness
h_ballast=0.6*hb_rel; % (KW51 - 04A) - GTN-D0-L201-03-006A

% 1) x-coordinates
% a) from left abutment to middle node
xn=[0,658,2313/3,658+401,2313/3*2,2313]; % [mm] 1 - 2
xn=[xn,xn(end)+1878*[1/2,1]]; % 2 - 3 
xn=[xn,xn(end)+3534*[1/2,1]]; % 3 - 4
xn=[xn,xn(end)+3533*[1/2,1]]; % 4 - 5 - connection at first diagonal
xn=[xn,xn(end)+3833*[1/2,1]]; % 5 - 6
xn=[xn,xn(end)+[3834/2,5900-3833,3834]]; % 6 - 7; additional node where 
                                         % diagonal elements connect to 
                                         % neutral axis of the arch
xn=[xn,xn(end)+3833*[1/2,1]]; % 7 - 8 - connection to second diagonal
for ind1=1:3 % 8 - 17    
    xn=[xn,xn(end)+3933*[1/2,1]];
    xn=[xn,xn(end)+3934*[1/2,1]];
    xn=[xn,xn(end)+3933*[1/2,1]]; % connection to diagonal 
end
clear ind1
xn=[xn,115000/2+658-13767-23400]; % 7 - 8 additional node where transversal 
                                  % stiffener connects to the arch
xn=[xn,115000/2+658-13767];       % 7 - 8 additional node where transversal 
                                  % stiffener connects to the arch
xn=sort(xn);

% b) from middle node to right abutment
xn=[xn,xn(end)+(xn(end)-xn(end-1:-1:1))];
xn=[xn,115000/2+658+3933/2+300];  % additional node where strain gauge
                                  % sgBD1718C is connected                       
xn=sort(xn).'/1000;  % ## [mm] => [m]

% 2) y-coordinates
yn=ones(size(xn));

% 3) Assemble nodes matrix
Nodes=[...
    ... % ## Main girder A (y = 0) - neutral axis (101XXX), XXX=000->REF
    101000,0,0,20;... 
    101000+(1:length(xn)).',xn,0*yn,0.2719*ones(size(xn));...
    ... % ## Main girder C - neutral axis (102XXX)
    102000,0,11.750,20;... 
    102000+(1:length(xn)).',xn,11.750*yn,0.2719*ones(size(xn));...
    ... % ## Arch A - Neutral axis 
    % (z = -0.0056679878x^2 + 0.659277672x + 0.341695) (103XXX)
    103000+(2:length(xn)-1).',xn(2:end-1),0*yn(2:end-1),...
        -0.0056679878*xn(2:end-1).^2+0.659277672*xn(2:end-1)+0.341695;  
    ... % ## Arch C - Neutral axis 
    % (z = -0.0056679878x^2 + 0.659277672x + 0.341695) (104XXX)
    104000+(2:length(xn)-1).',xn(2:end-1),11.750*yn(2:end-1),...
        -0.0056679878*xn(2:end-1).^2+0.659277672*xn(2:end-1)+0.341695]; 

% ## Nodes where diagonals connect A (105XXX)

% x-coordinate lower nodes where diagonals connect (on main girders) 
xdl=658+10600;
xdl=[xdl,xdl(end)+11500];
for ind1=1:6, xdl=[xdl,xdl(end)+11800]; end
clear ind1
xdl=[xdl,xdl(end)+11500];
xdl=xdl.'/1000;

% x-coordinate upper nodes where diagonals connect (on arch) 
xdu=658+10600+5900;
xdu=[xdu,xdu(end)+5600+5900];
for ind1=1:5, xdu=[xdu,xdu(end)+11800]; end
clear ind1
xdu=[xdu,xdu(end)+5600+5900];
xdu=xdu.'/1000;

% Node numbers lower nodes
ndl=zeros(size(xdl));
for ind1=1:length(ndl)
    iSel=(floor(Nodes(:,1)/1000)==101)&(Nodes(:,2)==xdl(ind1));
    ndl(ind1)=105000+Nodes(iSel,1)-floor(Nodes(iSel,1)/1000)*1000;
end
clear iSel ind1

% Node numbers upper nodes
ndu=zeros(size(xdu));
for ind1=1:length(ndu)
    iSel=(floor(Nodes(:,1)/1000)==103)&(Nodes(:,2)==xdu(ind1));
    ndu(ind1)=Nodes(iSel,1);
end
clear iSel ind1
Nodes=[Nodes;
    ndl,xdl,0*xdl,(80+693)/1000*ones(size(xdl))];

% ## Nodes where diagonals connect C (106XXX)
Nodes=[Nodes;
    ndl+1000,xdl,11.750*ones(size(xdl)),(80+693)/1000*ones(size(xdl))];

% ## Braces (107XXX - 110XXX)
Nebr=6; % Number of elements per brace
xbr=0.658+115/2+[-23400-13767,-13767,+13767,+13767+23400]/1000;
xbr=repmat(xbr,Nebr+1,1);
nbr=zeros(Nebr+1,4);
ybr=repmat(linspace(0,11.750,Nebr).',1,4);

% ## Node required for damage simulation
ybr=[ybr(1:end-1,:);ybr(end,:)-0.85;ybr(end,:)]; 
zbr=zeros(Nebr+1,4);
for ind1=1:4
    iSel=(floor(Nodes(:,1)/1000)==103)&...
        (round(Nodes(:,2)*10000)==round(xbr(1,ind1)*10000));
    nbr(1,ind1)=Nodes(iSel,1);
    iSel=(floor(Nodes(:,1)/1000)==104)&...
        (round(Nodes(:,2)*10000)==round(xbr(1,ind1)*10000));
    nbr(end,ind1)=Nodes(iSel,1);
    zbr(:,ind1)=Nodes(iSel,4);
    clear iSel
    nbr(2:end-1,ind1)=107000+(ind1-1)*1000+(1:Nebr-1).';
end
clear ind1

% ## Reference nodes for braces
xbr=[xbr;xbr(end,:)];
ybr=[ybr;ybr(1,:)];
zbr=[zbr;20*ones(1,4)];
nbr=[nbr;107000:1000:110000];
Nodes=[Nodes;
    reshape(nbr([2:end-2,end],:),numel(nbr([2:end-2,end],:)),1),...
    reshape(xbr([2:end-2,end],:),numel(xbr([2:end-2,end],:)),1),...
    reshape(ybr([2:end-2,end],:),numel(ybr([2:end-2,end],:)),1),...
    reshape(zbr([2:end-2,end],:),numel(zbr([2:end-2,end],:)),1)];

% ## Nodes on diagonals A at 0.8 m above the bottom end of the 
%    rectangular tube 
% (111XXX - Leaning towards Leuven)
% (112XXX - Leaning towards Brussels)
ndlb=[ndl(1:end-1)+6000;ndl(2:end)+7000];
xdlb=zeros(size(ndlb));
zdlb=zeros(size(ndlb));
zdlb_ref=(80+693)/1000;
xdlb(1)=xdl(1)+(2.383+0.8)/10.9394*5.900; % KW51-05A - node 9
zdlb(1)=zdlb_ref+(2.383+0.8)/10.9394*9.212;
xdlb(2)=xdl(2)+(2.252+0.8)/15.0151*5.900; % KW51-05A - node 8
zdlb(2)=zdlb_ref+(2.252+0.8)/15.0151*13.8074;
xdlb(3)=xdl(3)+(2.302+0.8)/17.9606*5.900; % KW51-05A - node 7
zdlb(3)=zdlb_ref+(2.302+0.8)/17.9606*16.9639;
xdlb(4)=xdl(4)+(2.404+0.8)/19.4585*5.900; % KW51-05A - node 6
zdlb(4)=zdlb_ref+(2.404+0.8)/19.4585*18.5425;
xdlb(5)=xdl(5)+(2.404+0.8)/19.4585*5.900; % KW51-05A - node 5
zdlb(5)=zdlb_ref+(2.404+0.8)/19.4585*18.5425;
xdlb(6)=xdl(6)+(2.395+0.8)/17.9606*5.900; % KW51-05A - node 6
zdlb(6)=zdlb_ref+(2.395+0.8)/17.9606*16.9639;
xdlb(7)=xdl(7)+(2.302+0.8)/15.0151*5.900; % KW51-05A - node 7
zdlb(7)=zdlb_ref+(2.302+0.8)/15.0151*13.8074;
xdlb(8)=xdl(8)+(2.399+0.8)/10.7806*5.600; % KW51-05A - node 8
zdlb(8)=zdlb_ref+(2.399+0.8)/10.7806*9.212;
xdlb(9)=xdl(2)-(2.399+0.8)/10.7806*5.600; % KW51-05A - node 8
zdlb(9)=zdlb_ref+(2.399+0.8)/10.7806*9.212;
xdlb(10)=xdl(3)-(2.302+0.8)/15.0151*5.900; % KW51-05A - node 7
zdlb(10)=zdlb_ref+(2.302+0.8)/15.0151*13.8074;
xdlb(11)=xdl(4)-(2.395+0.8)/17.9606*5.900; % KW51-05A - node 6
zdlb(11)=zdlb_ref+(2.395+0.8)/17.9606*16.9639;
xdlb(12)=xdl(5)-(2.404+0.8)/19.4585*5.900; % KW51-05A - node 5
zdlb(12)=zdlb_ref+(2.404+0.8)/19.4585*18.5425;
xdlb(13)=xdl(6)-(2.404+0.8)/19.4585*5.900; % KW51-05A - node 6
zdlb(13)=zdlb_ref+(2.404+0.8)/19.4585*18.5425;
xdlb(14)=xdl(7)-(2.302+0.8)/17.9606*5.900; % KW51-05A - node 7
zdlb(14)=zdlb_ref+(2.302+0.8)/17.9606*16.9639;
xdlb(15)=xdl(8)-(2.252+0.8)/15.0151*5.900; % KW51-05A - node 8
zdlb(15)=zdlb_ref+(2.252+0.8)/15.0151*13.8074;
xdlb(16)=xdl(9)-(2.383+0.8)/10.9394*5.900; % KW51-05A - node 9
zdlb(16)=zdlb_ref+(2.383+0.8)/10.9394*9.212;

% ## Reference nodes for diagonals
xdlb=[xdlb;0];
zdlb=[zdlb;20];
ndlb=[ndlb;111000];
Nodes=[Nodes;
    ndlb,xdlb,0*xdlb,zdlb];

% ## Nodes on diagonals A at 0.8 m below the top end of the 
%    rectangular tube 
% (113XXX - Leaning towards Leuven)
% (114XXX - Leaning towards Brussels)
ndub=[ndu+10000;ndu+11000];
xdub=zeros(size(ndub));
zdub=zeros(size(ndub));
zdub_ref=-0.0056679878*xdu.^2+0.659277672*xdu+0.341695;
zdub_ref=[zdub_ref;zdub_ref];
xdub(1)=xdu(1)-(3.484+0.8)/10.9394*5.900;    % KW51-04A - node 4
zdub(1)=zdub_ref(1)-(3.484+0.8)/10.9394*9.212;
xdub(2)=xdu(2)-(2.741+0.8)/15.0151*5.900;    % KW51-04A - node 3
zdub(2)=zdub_ref(2)-(2.741+0.8)/15.0151*13.8074;
xdub(3)=xdu(3)-(2.559+0.8)/17.9606*5.900;    % KW51-03 - node 2
zdub(3)=zdub_ref(3)-(2.559+0.8)/17.9606*16.9639;
xdub(4)=xdu(4)-(2.502+0.8)/19.4585*5.900;    % KW51-03 - node 1
zdub(4)=zdub_ref(4)-(2.502+0.8)/19.4585*18.5425;
xdub(5)=xdu(5)-(2.4136+0.8)/19.4585*5.900;   % KW51-03 - node 1
zdub(5)=zdub_ref(5)-(2.4136+0.8)/19.4585*18.5425;
xdub(6)=xdu(6)-(2.402+0.8)/17.9606*5.900;    % KW51-03 - node 2
zdub(6)=zdub_ref(6)-(2.402+0.8)/17.9606*16.9639;
xdub(7)=xdu(7)-(1.389+0.8)/15.0151*5.900;    % KW51-04A - node 3
zdub(7)=zdub_ref(7)-(1.389+0.8)/15.0151*13.8074;
xdub(8)=xdu(8)-(2.652+0.8)/10.7806*5.600;    % KW51-04A - node 4
zdub(8)=zdub_ref(8)-(2.652+0.8)/10.7806*9.212;
xdub(16)=xdu(8)+(3.484+0.8)/10.9394*5.900;   % KW51-04A - node 4
zdub(16)=zdub_ref(16)-(3.484+0.8)/10.9394*9.212;
xdub(15)=xdu(7)+(2.741+0.8)/15.0151*5.900;   % KW51-04A - node 3
zdub(15)=zdub_ref(15)-(2.741+0.8)/15.0151*13.8074;
xdub(14)=xdu(6)+(2.559+0.8)/17.9606*5.900;   % KW51-03 - node 2
zdub(14)=zdub_ref(14)-(2.559+0.8)/17.9606*16.9639;
xdub(13)=xdu(5)+(2.502+0.8)/19.4585*5.900;   % KW51-03 - node 1
zdub(13)=zdub_ref(13)-(2.502+0.8)/19.4585*18.5425;
xdub(12)=xdu(4)+(2.4136+0.8)/19.4585*5.900;  % KW51-03 - node 1
zdub(12)=zdub_ref(12)-(2.4136+0.8)/19.4585*18.5425;
xdub(11)=xdu(3)+(2.402+0.8)/17.9606*5.900;   % KW51-03 - node 2
zdub(11)=zdub_ref(11)-(2.402+0.8)/17.9606*16.9639;
xdub(10)=xdu(2)+(1.389+0.8)/15.0151*5.900;   % KW51-04A - node 3
zdub(10)=zdub_ref(10)-(1.389+0.8)/15.0151*13.8074;
xdub(9)=xdu(1)+(2.652+0.8)/10.7806*5.600;    % KW51-04A - node 4
zdub(9)=zdub_ref(9)-(2.652+0.8)/10.7806*9.212;

% ## Reference nodes for diagonals
xdub=[xdub;0];
zdub=[zdub;20];
ndub=[ndub;113000];
Nodes=[Nodes;
    ndub,xdub,0*xdub,zdub];

% ## Nodes on diagonals C at 0.8 m above the bottom end of the 
%    rectangular tube 
% (115XXX - Leaning towards Leuven)
% (116XXX - Leaning towards Brussels)
Nodes=[Nodes;
    ndlb+4000,xdlb,11.750*ones(size(xdlb)),zdlb];

% ## Nodes on diagonals C at 0.8 m below the top end of the 
%    rectangular tube 
% (117XXX - Leaning towards Leuven)
% (118XXX - Leaning towards Brussels)
Nodes=[Nodes;
    ndub+4000,xdub,11.750*ones(size(xdub)),zdub];

% ## U-shaped beams of the orthogonal bridge deck (121XXX-132XXX)
% Bottom assumed to coincide with the neutral axis of the bridge deck 
% z = 0.773 m + 0.238 m -> 2% inclination is disregarded (see KW51 -04A)
yu=(11750/2-375-5*750+(0:750:(11*750)))/1000;
for ind1=1:length(yu)
    Nodes=[Nodes;
        121000+ind1*1000,0,yu(ind1),20;... 
        121000+ind1*1000+(1:length(xn)).',xn,yu(ind1)*yn,1.011*ones(size(xn))];
end
clear ind1

% ## Transverse beams (141XXX - 173XXX)
xtr=[0,2313];
xtr=[xtr,xtr(end)+1878];
xtr=[xtr,xtr(end)+3534];
xtr=[xtr,xtr(end)+3533];
xtr=[xtr,xtr(end)+3833];
xtr=[xtr,xtr(end)+3834];
xtr=[xtr,xtr(end)+3833];
for ind1=1:3
    xtr=[xtr,xtr(end)+3933];
    xtr=[xtr,xtr(end)+3934];
    xtr=[xtr,xtr(end)+3933];
end
clear ind1
xtr=[xtr,xtr(end)+(xtr(end)-xtr(end-1:-1:1))];
xtr=xtr/1000;
xtr=repmat(xtr,14,1);
ytr=repmat([0,yu,11.75].',1,33);
ztr=0.2539*ones(size(xtr));
ntr=zeros(size(xtr));
for ind1=1:33
    for ind2=1:14
        ntr(ind2,ind1)=140000+ind1*1000+ind2;
    end
end
clear ind1 ind2

% ## Reference nodes for transverse beams
xtr=[xtr;xtr(end,:)];
ytr=[ytr;ytr(1,:)];
ztr=[ztr;20*ones(1,33)];
ntr=[ntr;floor(ntr(end,:)/1000)*1000];
Nodes=[Nodes;
    reshape(ntr,numel(ntr),1),...
    reshape(xtr,numel(xtr),1),...
    reshape(ytr,numel(ytr),1),...
    reshape(ztr,numel(ztr),1)];

% ## Nodes of the plate elements (181XXX - 196XXX)
% 2% slope of plate is neglected. 
% Plate level = top level U-shaped beams z = 0.773 m + 0.383 m
yp=[0,5.875-5.076,(11750/2-375-5*750+(0:750:(11*750)))/1000,5.875+5.076,11.750];
Nodes=[Nodes;
     180000+1000+(1:length(xn)).',xn,yp(1)*yn,1.156*ones(size(xn));
     180000+2000+(1:length(xn)).',xn,yp(2)*yn,1.156*ones(size(xn))];
for ind1=3:length(yp)-2
    Nodes=[Nodes;
        180000+ind1*1000+(1:length(xn)).',xn,[yp(ind1)-0.32/2]*yn,1.156*ones(size(xn));
        180500+ind1*1000+(1:length(xn)).',xn,[yp(ind1)+0.32/2]*yn,1.156*ones(size(xn))];
end
clear ind1
Nodes=[Nodes;
     180000+(length(yp)-1)*1000+(1:length(xn)).',xn,yp(end-1)*yn,1.156*ones(size(xn));
     180000+length(yp)*1000+(1:length(xn)).',xn,yp(end)*yn,1.156*ones(size(xn))];
 
% ## Nodes of the ballast layer (282XXX - 295XXX)
Nodes=[Nodes;
     280000+2000+(1:length(xn)).',xn,yp(2)*yn,(1.156+h_ballast/2)*ones(size(xn))];
for ind1=3:length(yp)-2
    Nodes=[Nodes;
        280000+ind1*1000+(1:length(xn)).',xn,[yp(ind1)-0.32/2]*yn,(1.156+h_ballast/2)*ones(size(xn));
        280500+ind1*1000+(1:length(xn)).',xn,[yp(ind1)+0.32/2]*yn,(1.156+h_ballast/2)*ones(size(xn))];
end
clear ind1
Nodes=[Nodes;
     280000+(length(yp)-1)*1000+(1:length(xn)).',xn,yp(end-1)*yn,(1.156+h_ballast/2)*ones(size(xn))];
 
% ## Nodes to apply boundary conditions (abutments)

% Hinge (constraint x, y, z) - Brussels, bottom (N1)
iSel1=(floor(Nodes(:,1)/1000)==101)&...
        (round(Nodes(:,2)*10000)==round(1.059*10000));
Nodes=[Nodes;1,Nodes(iSel1,2:3),-0.15];
clear iSel1 

% Constraint (constraint x, z) - Brussels, top (N2)
iSel2=(floor(Nodes(:,1)/1000)==102)&...
    (round(Nodes(:,2)*10000)==round(1.059*10000));
Nodes=[Nodes;2,Nodes(iSel2,2:3),-0.15];
clear iSel2 

% Constraint (constraint y, z) - Leuven, bottom (N3)
iSel3=(floor(Nodes(:,1)/1000)==101)&...
    (round(Nodes(:,2)*10000)==round((115000+658-401)/1000*10000));
Nodes=[Nodes;3,Nodes(iSel3,2:3),-0.15];
clear iSel3 

% Constraint (constraint z) - Leuven, top (N4)
iSel4=(floor(Nodes(:,1)/1000)==102)&...
    (round(Nodes(:,2)*10000)==round((115000+658-401)/1000*10000));
Nodes=[Nodes;4,Nodes(iSel4,2:3),-0.15];
clear iSel4 

% Clear obsolete variables
clear xbr xdk xdlb xdu xdub ybr yn ytr zbr zdlb zdlb_ref zdub zdub_ref ztr

% ## Plot nodes
if plot_elem==1 
    figure
    plotnodes(Nodes(Nodes(:,4)<20,:),'Numbering','off');
    % Plot reference nodes
    % hold on   
    % plotnodes(Nodes(Nodes(:,4)>=20,:),'Numbering','off','color','m','GCS','off');
end

%% TYPES DEFINITION

% Element types -> {EltTypID EltName}

Types={1 'beam';
       2 'truss';
       3 'shell4'};

%% SECTIONS DEFINITION

% Sections=[SecID A ky kz Ixx Iyy Izz]

% The following sections properties have been calculated from the geometry
% in the blueprints using ANSYS

Sections=...
    [1,7.688e-2,2.93,1.91,1.010e-4,1.440e-3,1.010e-2;... % Main girder 
     2,2.097e-1,2.04,2.84,6.610e-2,3.790e-2,5.430e-2;... % Arch 
     3,2.089e-2,2.35,2.30,5.960e-4,3.840e-4,3.760e-4;... % Diagonals 
     4,3.289e-2,2.00,2.00,5.770e-3,2.880e-3,2.880e-3;... % Braces 
     5,1.339e-2,2.21,2.68,3.060e-4,1.640e-4,2.900e-4;... % U-shaped beams 
     6,5.244e-2,2.65,2.08,2.760e-5,9.010e-4,5.630e-3;... % Transverse beams 
     7,0.015,NaN,NaN,NaN,NaN,NaN;...                     % Plate thickness 
     8,h_ballast,NaN,NaN,NaN,NaN,NaN;...                 % Ballast thickness
     9,2.089e-2*Ad_rel,2.35,2.30,5.960e-4,3.840e-4*Id_rel,3.760e-4;...  % Diagonals - connection
     10,2.089e-2*0.94,2.35,2.30,5.960e-4,3.840e-4*0.085,3.760e-4];      % Diagonals - connection (after retrofit)
 
%% MATERIALS DEFINITION

% Materials=[MatID E nu rho];
Materials=...
    [1,200e9*Es_rel,0.3,7850;... % Steel beam members
     2,200e9*Es_rel,0.3,7850;... % Steel plates
     3,550e6,0.3,1700;...        % Porphyry ballast layer 
     4,200e9*Es_rel,0.3,7850];   % Connection of diagonals with deck and arches

%% ELEMENTS DEFINITION

% Elements=[EltID TypID SecID MatID n1 n2 n3 n4]

Elements=[];

% ## Main girder A
for ind1=1:length(xn)-1
    Elements=[Elements;...
        101000+ind1,1,1,1,101000+ind1,101001+ind1,101000];
end

% ## Main girder C
for ind1=1:length(xn)-1
    Elements=[Elements;...
        102000+ind1,1,1,1,102000+ind1,102001+ind1,102000];
end

% ## Arch A
for ind1=2:length(xn)-2
    Elements=[Elements;...
        103000+ind1,1,2,1,103000+ind1,103001+ind1,101000];
end

% ## Arch C
for ind1=2:length(xn)-2
    Elements=[Elements;...
        104000+ind1,1,2,1,104000+ind1,104001+ind1,102000];
end

% ## Diagonals A (side of track B)
nd=ndl([1:8,2:9]);
nd=[nd,ndlb(1:end-1),ndub(1:end-1),ndu([1:8,1:8])];
nd=nd([16,8,15,7,14,6,13,5,12,4,11,3,10,2,9,1],:); 
ndref=ndlb(end);
for ind1=1:16

    if sum((floor(diags_box/100)==10)+((diags_box-floor(diags_box/100)*100)==ind1)==2)>0
        Elements=[Elements;... % Bottom
            105000+(ind1-1)*3+1,1,10,4,nd(ind1,1),nd(ind1,2),ndref];
    else
        Elements=[Elements;... % Bottom
            105000+(ind1-1)*3+1,1,9,4,nd(ind1,1),nd(ind1,2),ndref];
    end
    
    Elements=[Elements;... % Middle
        105000+(ind1-1)*3+2,1,3,1,nd(ind1,2),nd(ind1,3),ndref];
    
    if sum((floor(diags_box/100)==20)+((diags_box-floor(diags_box/100)*100)==ind1)==2)>0
        Elements=[Elements;... % Top
            105000+(ind1-1)*3+3,1,10,4,nd(ind1,3),nd(ind1,4),ndref];
    else
        Elements=[Elements;... % Top
            105000+(ind1-1)*3+3,1,9,4,nd(ind1,3),nd(ind1,4),ndref];
    end

end

% ## Diagonals C (side of track A)
for ind1=1:16

    if sum((floor(diags_box/100)==30)+((diags_box-floor(diags_box/100)*100)==ind1)==2)>0
        Elements=[Elements;... % Bottom
            106000+(ind1-1)*3+1,1,10,4,nd(ind1,1)+1000,nd(ind1,2)+4000,ndref+4000];
    else
        Elements=[Elements;... % Bottom
            106000+(ind1-1)*3+1,1,9,4,nd(ind1,1)+1000,nd(ind1,2)+4000,ndref+4000];        
    end

    Elements=[Elements;... % Middle
        106000+(ind1-1)*3+2,1,3,1,nd(ind1,2)+4000,nd(ind1,3)+4000,ndref+4000];
    
    if sum((floor(diags_box/100)==40)+((diags_box-floor(diags_box/100)*100)==ind1)==2)>0
        Elements=[Elements;... % Top
            106000+(ind1-1)*3+3,1,10,4,nd(ind1,3)+4000,nd(ind1,4)+1000,ndref+4000];
    else
        Elements=[Elements;... % Top
            106000+(ind1-1)*3+3,1,9,4,nd(ind1,3)+4000,nd(ind1,4)+1000,ndref+4000];
    end

end

% ## Braces
for ind1=1:4
    for ind2=1:Nebr
        Elements=[Elements;...
            107000+(ind1-1)*1000+ind2,1,4,1,...
            nbr(ind2,ind1),nbr(ind2+1,ind1),nbr(end,ind1)];
    end
end

% ## U-shaped beams of the orthogonal bridge deck 
for ind1=1:12
    for ind2=1:length(xn)-1
        Elements=[Elements;...
            121000+ind1*1000+ind2,1,5,1,...
            121000+ind1*1000+[ind2,ind2+1,0]];
    end
end

% ## Transverse beams
for ind1=1:33
    for ind2=1:13
        Elements=[Elements;...
            140000+ind2*1000+ind1,1,6,1,...
            ntr(ind2,ind1),ntr(ind2+1,ind1),ntr(end,ind1)];
    end
end

% ## Plates (bridge deck)
Elements=[Elements,NaN*ones(size(Elements,1),1)];
for ind2=1:length(xn)-1
    Elements=[Elements;...
        180000+1000+ind2,3,7,2,...
            180000+1000+[ind2,ind2+1],180000+2000+[ind2+1,ind2]];
    Elements=[Elements;...
        180000+2000+ind2,3,7,2,...
            180000+2000+[ind2,ind2+1],180000+3000+[ind2+1,ind2]];
end
for ind1=3:length(yp)-2
    for ind2=1:length(xn)-1
        Elements=[Elements;...
            180000+ind1*1000+ind2,3,7,2,180500+ind1*1000+[ind2,ind2+1],...
            180000+(ind1+1)*1000+[ind2+1,ind2]];
    end
end
for ind2=1:length(xn)-1
    Elements=[Elements;...
        180000+(length(yp)-1)*1000+ind2,3,7,2,...
            180000+(length(yp)-1)*1000+[ind2,ind2+1],...
            180000+(length(yp))*1000+[ind2+1,ind2]];
end

% ## Ballast
for ind2=1:length(xn)-1
    Elements=[Elements;...
        280000+2000+ind2,3,8,3,...
            280000+2000+[ind2,ind2+1],...
            280000+3000+[ind2+1,ind2]];
end
for ind1=3:length(yp)-2
    for ind2=1:length(xn)-1
        Elements=[Elements;...
            280000+ind1*1000-500+ind2,3,8,3,...
            280000+(ind1)*1000+[ind2+1,ind2],...
            280500+(ind1)*1000+[ind2,ind2+1]];
        Elements=[Elements;...
            280000+ind1*1000+ind2,3,8,3,...
            280500+ind1*1000+[ind2,ind2+1],...
            280000+(ind1+1)*1000+[ind2+1,ind2]];
    end
end

% ## Plot elements
if plot_elem==1 
    hold on
    plotelem(Nodes,Elements(Elements(:,4)==1,:),Types,'Numbering','off','GCS','off');
    plotelem(Nodes,Elements(Elements(:,4)==2,:),Types,'Numbering','off','color','r','GCS','off');
    plotelem(Nodes,Elements(Elements(:,4)==3,:),Types,'Numbering','off','color','c','GCS','off');
    plotelem(Nodes,Elements(Elements(:,4)==4,:),Types,'Numbering','off','color','b','GCS','off');
end

%% ASSEMBLY OF STIFFNESS AND MASS MATRIX

% Degrees of freedom vector
DOF=getdof(Elements,Types);

% Add nodes to apply boundary conditions
DOF=[DOF;(1:4).'+0.01;(1:4).'+0.02;(1:4).'+0.03];

% Assemble the stiffness and mass matrix
[K,M]=asmkm(Nodes,Elements,Types,Sections,Materials,DOF);

%% ACCOUNT FOR THE MASS OF THE STEEL BOXES ADDED DURING THE RETROFIT

mbox=333; % mass boxes = 333 kg

if addmass==1&&isempty(diags_box) 

    % ## Add mass for all diagonals, top and bottom
    
    % Track B, bottom
    for ind1=1:size(nd,1)
        seldof=[nd(ind1,2)+0.01;nd(ind1,2)+0.02;nd(ind1,2)+0.03];
        [~,I]=selectdof(DOF,seldof);
        for ind2=1:length(I)
            M(I(ind2),I(ind2))=M(I(ind2),I(ind2))+mbox;
        end
    end
    
    % Track B, top
    for ind1=1:size(nd,1)
        seldof=[nd(ind1,3)+0.01;nd(ind1,3)+0.02;nd(ind1,3)+0.03];
        [~,I]=selectdof(DOF,seldof);
        for ind2=1:length(I)
            M(I(ind2),I(ind2))=M(I(ind2),I(ind2))+mbox;
        end
    end
    
    % Track A, bottom
    for ind1=1:size(nd,1)
        seldof=[nd(ind1,2)+4000+0.01;nd(ind1,2)+4000+0.02;nd(ind1,2)+4000+0.03];
        [~,I]=selectdof(DOF,seldof);
        for ind2=1:length(I)
            M(I(ind2),I(ind2))=M(I(ind2),I(ind2))+mbox;
        end
    end
    
    % Track A, top
    for ind1=1:size(nd,1)
        seldof=[nd(ind1,3)+4000+0.01;nd(ind1,3)+4000+0.02;nd(ind1,3)+4000+0.03];
        [~,I]=selectdof(DOF,seldof);
        for ind2=1:length(I)
            M(I(ind2),I(ind2))=M(I(ind2),I(ind2))+mbox;
        end
    end
      
    clear I seldof Nodes_boxes
    
elseif isempty(diags_box)==0

    % ## Add mass for a selection of diagonals (stiffness will be accounted
    % for in a later stage)
    
    % Track B, bottom
    for ind1=1:size(nd,1)
        if sum((floor(diags_box/100)==10)+((diags_box-floor(diags_box/100)*100)==ind1)==2)>0
            seldof=[nd(ind1,2)+0.01;nd(ind1,2)+0.02;nd(ind1,2)+0.03];
            [~,I]=selectdof(DOF,seldof);
            for ind2=1:length(I)
                M(I(ind2),I(ind2))=M(I(ind2),I(ind2))+mbox;
            end
        end
    end
    
    % Track B, top
    for ind1=1:size(nd,1)
        if sum((floor(diags_box/100)==20)+((diags_box-floor(diags_box/100)*100)==ind1)==2)>0
            seldof=[nd(ind1,3)+0.01;nd(ind1,3)+0.02;nd(ind1,3)+0.03];
            [~,I]=selectdof(DOF,seldof);
            for ind2=1:length(I)
                M(I(ind2),I(ind2))=M(I(ind2),I(ind2))+mbox;
            end
        end
    end
    
    % Track A, bottom
    for ind1=1:size(nd,1)
        if sum((floor(diags_box/100)==30)+((diags_box-floor(diags_box/100)*100)==ind1)==2)>0
            seldof=[nd(ind1,2)+4000+0.01;nd(ind1,2)+4000+0.02;nd(ind1,2)+4000+0.03];
            [~,I]=selectdof(DOF,seldof);
            for ind2=1:length(I)
                M(I(ind2),I(ind2))=M(I(ind2),I(ind2))+mbox;
            end
        end
    end
    
    % Track A, top
    for ind1=1:size(nd,1)
        if sum((floor(diags_box/100)==40)+((diags_box-floor(diags_box/100)*100)==ind1)==2)>0
            seldof=[nd(ind1,3)+4000+0.01;nd(ind1,3)+4000+0.02;nd(ind1,3)+4000+0.03];
            [~,I]=selectdof(DOF,seldof);
            for ind2=1:length(I)
                M(I(ind2),I(ind2))=M(I(ind2),I(ind2))+mbox;
            end
        end
    end
      
    clear I seldof Nodes_boxes
    
end

%% ADD CONSTRAINT EQUATIONS

% Help from the tconstr function 
% [T,Q0,MasterDOF]=tconstr(Constr,DOF)
%  returns matrices to apply constraint equations to the stiffness and mass 
%  matrix and the load vector: Kr=T.'*K*T, Mr=T.'*M*T and Fr=T.'*(F-K*Q0). 
%  The original displacement vector is computed using U=T*Ur+Q0. 
%  Constr     Constraint equation
%      Constant=CoefS*SlaveDOF+CoefM1*MasterDOF1+CoefM2*MasterDOF2+...
%     [Constant CoefS SlaveDOF CoefM1 MasterDOF1 CoefM2 MasterDOF2 ...]

% Overview of the steps where the constraints are applied
%       Master          Slave
% 1.1	Plate           Ballast
% ==> T1,MasterDOF1
% 2.1	Main girders	Diagonals
% 2.2	Main girders	Arch begin and end
% 2.3	Main girders	Plate
% 2.4	U-shaped beams	Plate
% 2.5	U-shaped beams	Transverse beams
% ==> T2,MasterDOF2
% 3.1   Slave main girders to bearings 
% ==> T3,MasterDOF3
% 4.1   Apply boundary conditions
% ==> T4,MasterDOF4
% T=T1*T2*T3*T4;
    
% ## Slave ballast (slave) to plates (master) 
Constr=[];
for ind1=2:(length(yp)-1)
    for ind2=1:length(xn)

        % Find node on plates (master node)
        iSel1=(floor(Nodes(:,1)/1000)==180+ind1)&...
            (Nodes(:,4)<20)&...
            (round(Nodes(:,2)*10000)==round(xn(ind2)*10000));
        Node1=Nodes(iSel1,:);   

        % Find node on ballast (slave nodes)
        iSel2=(floor(Nodes(:,1)/1000)==280+ind1)&...
            (Nodes(:,4)<20)&...
            (round(Nodes(:,2)*10000)==round(xn(ind2)*10000));
        Node2=Nodes(iSel2,:); 
        
        % Rotations node 2.1
        Constr=[Constr;
                0,-1,Node2(1,1)+0.04,1,Node1(1,1)+0.04,NaN,NaN,NaN,NaN;
                0,-1,Node2(1,1)+0.05,1,Node1(1,1)+0.05,NaN,NaN,NaN,NaN;
                0,-1,Node2(1,1)+0.06,1,Node1(1,1)+0.06,NaN,NaN,NaN,NaN];

        % Translations node 2.1
        dx=(Node2(1,2)-Node1(1,2));
        dy=(Node2(1,3)-Node1(1,3)); 
        dz=(Node2(1,4)-Node1(1,4));
        
        Constr=[Constr;
                0,-1,Node2(1,1)+0.01,1,Node1(1,1)+0.01,+dz,Node1(1,1)+0.05,NaN,NaN;
                0,-1,Node2(1,1)+0.02,1,Node1(1,1)+0.02,-dz,Node1(1,1)+0.04,NaN,NaN;
                0,-1,Node2(1,1)+0.03,1,Node1(1,1)+0.03,NaN,NaN,NaN,NaN];
        
        if ind1>2&&ind1<(length(yp)-1)

            % Rotations node 2.2
            Constr=[Constr;
                    0,-1,Node2(2,1)+0.04,1,Node1(2,1)+0.04,NaN,NaN,NaN,NaN;
                    0,-1,Node2(2,1)+0.05,1,Node1(2,1)+0.05,NaN,NaN,NaN,NaN;
                    0,-1,Node2(2,1)+0.06,1,Node1(2,1)+0.06,NaN,NaN,NaN,NaN];

            % Translations node 2.2
            dx=(Node2(2,2)-Node1(2,2));
            dy=(Node2(2,3)-Node1(2,3)); 
            dz=(Node2(2,4)-Node1(2,4));
            
            Constr=[Constr;
                    0,-1,Node2(2,1)+0.01,1,Node1(2,1)+0.01,+dz,Node1(2,1)+0.05,NaN,NaN;
                    0,-1,Node2(2,1)+0.02,1,Node1(2,1)+0.02,-dz,Node1(2,1)+0.04,NaN,NaN;
                    0,-1,Node2(2,1)+0.03,1,Node1(2,1)+0.03,NaN,NaN,NaN,NaN];

        end
    end
end

[T1,~,MasterDOF1]=tconstr(Constr,DOF);
clear Constr

% ## Slave diagonals (slave) to main girders (master)
Constr=[];
for ind1=1:length(xdl)
    
    % Main girder A
    iSel1=(floor(Nodes(:,1)/1000)==101)&...
        (round(Nodes(:,2)*10000)==round(xdl(ind1)*10000));
    Node1=Nodes(iSel1,:);   % Node on main girder (master node)
    iSel2=(floor(Nodes(:,1)/1000)==105)&...
        (round(Nodes(:,2)*10000)==round(xdl(ind1)*10000));
    Node2=Nodes(iSel2,:);   % Node on diagonals (slave node)
    dx=(Node2(2)-Node1(2));
    dy=(Node2(3)-Node1(3)); 
    dz=(Node2(4)-Node1(4));

    Constr=[Constr;
            0,-1,Node2(1)+0.01,1,Node1(1)+0.01,+dz,Node1(1)+0.05,NaN,NaN;
            0,-1,Node2(1)+0.02,1,Node1(1)+0.02,-dz,Node1(1)+0.04,NaN,NaN;
            0,-1,Node2(1)+0.03,1,Node1(1)+0.03,NaN,NaN,NaN,NaN;
            0,-1,Node2(1)+0.04,1,Node1(1)+0.04,NaN,NaN,NaN,NaN;
            0,-1,Node2(1)+0.05,1,Node1(1)+0.05,NaN,NaN,NaN,NaN;
            0,-1,Node2(1)+0.06,1,Node1(1)+0.06,NaN,NaN,NaN,NaN];
    
    % Main girder C
    iSel1=(floor(Nodes(:,1)/1000)==102)&...
        (round(Nodes(:,2)*10000)==round(xdl(ind1)*10000));
    Node1=Nodes(iSel1,:);   % Node on main girder (master node)
    iSel2=(floor(Nodes(:,1)/1000)==106)&...
        (round(Nodes(:,2)*10000)==round(xdl(ind1)*10000));
    Node2=Nodes(iSel2,:);   % Node on diagonals (slave node)
    dx=(Node2(2)-Node1(2));
    dy=(Node2(3)-Node1(3)); 
    dz=(Node2(4)-Node1(4));

    Constr=[Constr;
            0,-1,Node2(1)+0.01,1,Node1(1)+0.01,+dz,Node1(1)+0.05,NaN,NaN;
            0,-1,Node2(1)+0.02,1,Node1(1)+0.02,-dz,Node1(1)+0.04,NaN,NaN;
            0,-1,Node2(1)+0.03,1,Node1(1)+0.03,NaN,NaN,NaN,NaN;
            0,-1,Node2(1)+0.04,1,Node1(1)+0.04,NaN,NaN,NaN,NaN;
            0,-1,Node2(1)+0.05,1,Node1(1)+0.05,NaN,NaN,NaN,NaN;
            0,-1,Node2(1)+0.06,1,Node1(1)+0.06,NaN,NaN,NaN,NaN];
    
end

% ## Slave arch begin and end (slave) to main girders (master)

% Main girder A, begin
iSel1=(floor(Nodes(:,1)/1000)==101)&...
    (round(Nodes(:,2)*10000)==round(xn(2)*10000));
Node1=Nodes(iSel1,:);   % Node on main girder (master node)
iSel2=(floor(Nodes(:,1)/1000)==103)&...
    (round(Nodes(:,2)*10000)==round(xn(2)*10000));
Node2=Nodes(iSel2,:);   % Node on arch (slave node)
dx=(Node2(2)-Node1(2));
dy=(Node2(3)-Node1(3)); 
dz=(Node2(4)-Node1(4));

Constr=[Constr;
        0,-1,Node2(1)+0.01,1,Node1(1)+0.01,+dz,Node1(1)+0.05,NaN,NaN;
        0,-1,Node2(1)+0.02,1,Node1(1)+0.02,-dz,Node1(1)+0.04,NaN,NaN;
        0,-1,Node2(1)+0.03,1,Node1(1)+0.03,NaN,NaN,NaN,NaN;
        0,-1,Node2(1)+0.04,1,Node1(1)+0.04,NaN,NaN,NaN,NaN;
        0,-1,Node2(1)+0.05,1,Node1(1)+0.05,NaN,NaN,NaN,NaN;
        0,-1,Node2(1)+0.06,1,Node1(1)+0.06,NaN,NaN,NaN,NaN];

% Main girder A, end
iSel1=(floor(Nodes(:,1)/1000)==101)&...
    (round(Nodes(:,2)*10000)==round(xn(end-1)*10000));
Node1=Nodes(iSel1,:);   % Node on main girder (master node)
iSel2=(floor(Nodes(:,1)/1000)==103)&...
    (round(Nodes(:,2)*10000)==round(xn(end-1)*10000));
Node2=Nodes(iSel2,:);   % Node on arch (slave node)
dx=(Node2(2)-Node1(2));
dy=(Node2(3)-Node1(3)); 
dz=(Node2(4)-Node1(4));

Constr=[Constr;
        0,-1,Node2(1)+0.01,1,Node1(1)+0.01,+dz,Node1(1)+0.05,NaN,NaN;
        0,-1,Node2(1)+0.02,1,Node1(1)+0.02,-dz,Node1(1)+0.04,NaN,NaN;
        0,-1,Node2(1)+0.03,1,Node1(1)+0.03,NaN,NaN,NaN,NaN;
        0,-1,Node2(1)+0.04,1,Node1(1)+0.04,NaN,NaN,NaN,NaN;
        0,-1,Node2(1)+0.05,1,Node1(1)+0.05,NaN,NaN,NaN,NaN;
        0,-1,Node2(1)+0.06,1,Node1(1)+0.06,NaN,NaN,NaN,NaN];

% Main girder C, begin
iSel1=(floor(Nodes(:,1)/1000)==102)&...
    (round(Nodes(:,2)*10000)==round(xn(2)*10000));
Node1=Nodes(iSel1,:);   % Node on main girder (master node)
iSel2=(floor(Nodes(:,1)/1000)==104)&...
    (round(Nodes(:,2)*10000)==round(xn(2)*10000));
Node2=Nodes(iSel2,:);   % Node on arch (slave node)
dx=(Node2(2)-Node1(2));
dy=(Node2(3)-Node1(3)); 
dz=(Node2(4)-Node1(4));

Constr=[Constr;
        0,-1,Node2(1)+0.01,1,Node1(1)+0.01,+dz,Node1(1)+0.05,NaN,NaN;
        0,-1,Node2(1)+0.02,1,Node1(1)+0.02,-dz,Node1(1)+0.04,NaN,NaN;
        0,-1,Node2(1)+0.03,1,Node1(1)+0.03,NaN,NaN,NaN,NaN;
        0,-1,Node2(1)+0.04,1,Node1(1)+0.04,NaN,NaN,NaN,NaN;
        0,-1,Node2(1)+0.05,1,Node1(1)+0.05,NaN,NaN,NaN,NaN;
        0,-1,Node2(1)+0.06,1,Node1(1)+0.06,NaN,NaN,NaN,NaN];

% Main girder C, end
iSel1=(floor(Nodes(:,1)/1000)==102)&...
    (round(Nodes(:,2)*10000)==round(xn(end-1)*10000));
Node1=Nodes(iSel1,:);   % Node on main girder (master node)
iSel2=(floor(Nodes(:,1)/1000)==104)&...
    (round(Nodes(:,2)*10000)==round(xn(end-1)*10000));
Node2=Nodes(iSel2,:);   % Node on arch (slave node)
dx=(Node2(2)-Node1(2));
dy=(Node2(3)-Node1(3)); 
dz=(Node2(4)-Node1(4));

Constr=[Constr;
        0,-1,Node2(1)+0.01,1,Node1(1)+0.01,+dz,Node1(1)+0.05,NaN,NaN;
        0,-1,Node2(1)+0.02,1,Node1(1)+0.02,-dz,Node1(1)+0.04,NaN,NaN;
        0,-1,Node2(1)+0.03,1,Node1(1)+0.03,NaN,NaN,NaN,NaN;
        0,-1,Node2(1)+0.04,1,Node1(1)+0.04,NaN,NaN,NaN,NaN;
        0,-1,Node2(1)+0.05,1,Node1(1)+0.05,NaN,NaN,NaN,NaN;
        0,-1,Node2(1)+0.06,1,Node1(1)+0.06,NaN,NaN,NaN,NaN];
    
% ## Slave plates (slave) to main girder A (master)
for ind2=1:length(xn)
    
    % Find node on main girder (master node)
    iSel1=(floor(Nodes(:,1)/1000)==101)&...
        (Nodes(:,4)<20)&...
        (round(Nodes(:,2)*10000)==round(xn(ind2)*10000));
    Node1=Nodes(iSel1,:);   

    % Find node on plates (slave nodes)
    iSel2=(floor(Nodes(:,1)/1000)==181)&...
        (Nodes(:,4)<20)&...
        (round(Nodes(:,2)*10000)==round(xn(ind2)*10000));
    Node2=Nodes(iSel2,:); 

    % Rotations 
    Constr=[Constr;
            0,-1,Node2(1,1)+0.04,1,Node1(1,1)+0.04,NaN,NaN,NaN,NaN;
            0,-1,Node2(1,1)+0.05,1,Node1(1,1)+0.05,NaN,NaN,NaN,NaN;
            0,-1,Node2(1,1)+0.06,1,Node1(1,1)+0.06,NaN,NaN,NaN,NaN];

    % Translations node 2.2
    dx=(Node2(1,2)-Node1(2));
    dy=(Node2(1,3)-Node1(3)); 
    dz=(Node2(1,4)-Node1(4));
    
    Constr=[Constr;
            0,-1,Node2(1,1)+0.01,1,Node1(1,1)+0.01,+dz,Node1(1,1)+0.05,NaN,NaN;
            0,-1,Node2(1,1)+0.02,1,Node1(1,1)+0.02,-dz,Node1(1,1)+0.04,NaN,NaN;
            0,-1,Node2(1,1)+0.03,1,Node1(1,1)+0.03,NaN,NaN,NaN,NaN];
    
end

% ## Slave plates (slave) to main girder C (master)
for ind2=1:length(xn)
    
    % Find node on main girder (master node)
    iSel1=(floor(Nodes(:,1)/1000)==102)&...
        (Nodes(:,4)<20)&...
        (round(Nodes(:,2)*10000)==round(xn(ind2)*10000));
    Node1=Nodes(iSel1,:);   

    % Find node on plates (slave nodes)
    iSel2=(floor(Nodes(:,1)/1000)==196)&...
        (Nodes(:,4)<20)&...
        (round(Nodes(:,2)*10000)==round(xn(ind2)*10000));
    Node2=Nodes(iSel2,:); 

    % Rotations 
    Constr=[Constr;
            0,-1,Node2(1,1)+0.04,1,Node1(1,1)+0.04,NaN,NaN,NaN,NaN;
            0,-1,Node2(1,1)+0.05,1,Node1(1,1)+0.05,NaN,NaN,NaN,NaN;
            0,-1,Node2(1,1)+0.06,1,Node1(1,1)+0.06,NaN,NaN,NaN,NaN];

    % Translations 
    dx=(Node2(1,2)-Node1(2));
    dy=(Node2(1,3)-Node1(3)); 
    dz=(Node2(1,4)-Node1(4));

    Constr=[Constr;
            0,-1,Node2(1,1)+0.01,1,Node1(1,1)+0.01,+dz,Node1(1,1)+0.05,NaN,NaN;
            0,-1,Node2(1,1)+0.02,1,Node1(1,1)+0.02,-dz,Node1(1,1)+0.04,NaN,NaN;
            0,-1,Node2(1,1)+0.03,1,Node1(1,1)+0.03,NaN,NaN,NaN,NaN];
    
end

% ## Slave plates (slave) to U-shaped beams (master)
for ind1=1:length(yu)
    for ind2=1:length(xn)

        % Find node on U-shape beam (master node)
        iSel1=(floor(Nodes(:,1)/1000)==121+ind1)&...
            (Nodes(:,4)<20)&...
            (round(Nodes(:,2)*10000)==round(xn(ind2)*10000));
        Node1=Nodes(iSel1,:);   

        % Find node on plates (slave nodes)
        iSel2=(floor(Nodes(:,1)/1000)==181+ind1+1)&...
            (Nodes(:,4)<20)&...
            (round(Nodes(:,2)*10000)==round(xn(ind2)*10000));
        Node2=Nodes(iSel2,:); 
        
        % Rotations node 2.1
        Constr=[Constr;
                0,-1,Node2(1,1)+0.04,1,Node1(1,1)+0.04,NaN,NaN,NaN,NaN;
                0,-1,Node2(1,1)+0.05,1,Node1(1,1)+0.05,NaN,NaN,NaN,NaN;
                0,-1,Node2(1,1)+0.06,1,Node1(1,1)+0.06,NaN,NaN,NaN,NaN];

        % Translations node 2.1
        dx=(Node2(1,2)-Node1(2));
        dy=(Node2(1,3)-Node1(3)); 
        dz=(Node2(1,4)-Node1(4));

        Constr=[Constr;
                0,-1,Node2(1,1)+0.01,1,Node1(1,1)+0.01,+dz,Node1(1,1)+0.05,-dy,Node1(1,1)+0.06;
                0,-1,Node2(1,1)+0.02,1,Node1(1,1)+0.02,-dz,Node1(1,1)+0.04,NaN,NaN;
                0,-1,Node2(1,1)+0.03,1,Node1(1,1)+0.03,+dy,Node1(1,1)+0.04,NaN,NaN];
        
        % Rotations node 2.2
        Constr=[Constr;
                0,-1,Node2(2,1)+0.04,1,Node1(1,1)+0.04,NaN,NaN,NaN,NaN;
                0,-1,Node2(2,1)+0.05,1,Node1(1,1)+0.05,NaN,NaN,NaN,NaN;
                0,-1,Node2(2,1)+0.06,1,Node1(1,1)+0.06,NaN,NaN,NaN,NaN];

        % Translations node 2.2
        dx=(Node2(2,2)-Node1(2));
        dy=(Node2(2,3)-Node1(3)); 
        dz=(Node2(2,4)-Node1(4));

        Constr=[Constr;
                0,-1,Node2(2,1)+0.01,1,Node1(1,1)+0.01,+dz,Node1(1,1)+0.05,-dy,Node1(1,1)+0.06;
                0,-1,Node2(2,1)+0.02,1,Node1(1,1)+0.02,-dz,Node1(1,1)+0.04,NaN,NaN;
                0,-1,Node2(2,1)+0.03,1,Node1(1,1)+0.03,+dy,Node1(1,1)+0.04,NaN,NaN];

    end
end

% ## Slave transverse beams (slave) to main girders and U-shaped beams 
% (master)
for ind1=1:size(ntr,1)-1
    for ind2=1:size(ntr,2)
        
        if ind1==1 % Slave to main girder A
            
            % Find node on main girder A (master node)
            iSel1=(floor(Nodes(:,1)/1000)==101)&...
                (Nodes(:,4)<20)&...
                (round(Nodes(:,2)*10000)==round(xtr(ind1,ind2)*10000));
            
        elseif ind1==size(ntr,1)-1 % Slave to main girder C
            
            % Find node on main girder C (master node)
            iSel1=(floor(Nodes(:,1)/1000)==102)&...
                (Nodes(:,4)<20)&...
                (round(Nodes(:,2)*10000)==round(xtr(ind1,ind2)*10000));
            
        else % Slave to U-shaped beam
            
            % Find node on U-shaped beam (master node)
            iSel1=(floor(Nodes(:,1)/1000)==(120+ntr(ind1,ind2)-floor(ntr(ind1,ind2)/1000)*1000))&...
                (Nodes(:,4)<20)&...
                (round(Nodes(:,2)*10000)==round(xtr(ind1,ind2)*10000));

        end
         
        Node1=Nodes(iSel1,:);  

        % Find node on transverse beam (slave nodes)
        Node2=Nodes(Nodes(:,1)==ntr(ind1,ind2),:); 

        % Rotations 
        Constr=[Constr;
                0,-1,Node2(1,1)+0.04,1,Node1(1,1)+0.04,NaN,NaN,NaN,NaN;
                0,-1,Node2(1,1)+0.05,1,Node1(1,1)+0.05,NaN,NaN,NaN,NaN;
                0,-1,Node2(1,1)+0.06,1,Node1(1,1)+0.06,NaN,NaN,NaN,NaN];

        % Translations 
        dx=(Node2(1,2)-Node1(2));
        dy=(Node2(1,3)-Node1(3)); 
        dz=(Node2(1,4)-Node1(4));

        Constr=[Constr;
                0,-1,Node2(1,1)+0.01,1,Node1(1,1)+0.01,+dz,Node1(1,1)+0.05,NaN,NaN;
                0,-1,Node2(1,1)+0.02,1,Node1(1,1)+0.02,-dz,Node1(1,1)+0.04,NaN,NaN;
                0,-1,Node2(1,1)+0.03,1,Node1(1,1)+0.03,NaN,NaN,NaN,NaN];
        
    end
end

[T2,~,MasterDOF2]=tconstr(Constr,MasterDOF1);
clear Constr

% ## Slave main girders (slave) to bearings (master)
% For bearing numbers (N1 - N4), check 'overview.pdf'

% Hinge (constraint x, y, z) - Brussels, bottom
iSel1=(floor(Nodes(:,1)/1000)==101)&...
        (round(Nodes(:,2)*10000)==round(1.059*10000));
Node1=Nodes(iSel1,1);
dz_bc=Nodes(iSel1,4)-Nodes(Nodes(:,1)==1,4);
Constr1x=[0,1,Node1+0.01,-dz_bc,Node1+0.05,-1,1.01];
Constr1y=[0,1,Node1+0.02,+dz_bc,Node1+0.04,-1,1.02];
Constr1z=[0,1,Node1+0.03,-1,1.03];
Constr=[Constr1x;Constr1y;Constr1z,0,0];
seldof=1+[0.01,0.02,0.03];
clear iSel1 Node1

% Constraint (constraint x, z) - Brussels, top
iSel2=(floor(Nodes(:,1)/1000)==102)&...
    (round(Nodes(:,2)*10000)==round(1.059*10000));
Node2=Nodes(iSel2,1);
Constr2x=[0,1,Node2+0.01,-dz_bc,Node2+0.05,-1,2.01];
Constr2y=[0,1,Node2+0.02,+dz_bc,Node2+0.04,-1,2.02];
Constr2z=[0,1,Node2+0.03,-1,2.03];
Constr=[Constr;Constr2x;Constr2y;Constr2z,0,0];
seldof=[seldof,2+[0.01,0.03]];
clear iSel2 Node2

% Constraint (constraint y, z) - Leuven, bottom
iSel3=(floor(Nodes(:,1)/1000)==101)&...
    (round(Nodes(:,2)*10000)==round((115000+658-401)/1000*10000));
Node3=Nodes(iSel3,1);
Constr3x=[0,1,Node3+0.01,-dz_bc,Node3+0.05,-1,3.01];
Constr3y=[0,1,Node3+0.02,+dz_bc,Node3+0.04,-1,3.02];
Constr3z=[0,1,Node3+0.03,-1,3.03];
Constr=[Constr;Constr3x;Constr3y;Constr3z,0,0];
if fixed==0
    seldof=[seldof,3+[0.02,0.03]];
elseif fixed==1
    seldof=[seldof,3+[0.01,0.02,0.03]];
end
clear iSel3 

% Constraint (constraint z) - Leuven, top
iSel4=(floor(Nodes(:,1)/1000)==102)&...
    (round(Nodes(:,2)*10000)==round((115000+658-401)/1000*10000));
Node4=Nodes(iSel4,1);
Constr4x=[0,1,Node4+0.01,-dz_bc,Node4+0.05,-1,4.01];
Constr4y=[0,1,Node4+0.02,+dz_bc,Node4+0.04,-1,4.02];
Constr4z=[0,1,Node4+0.03,-1,4.03];
Constr=[Constr;Constr4x;Constr4y;Constr4z,0,0];
if fixed==0
    seldof=[seldof,4+0.03];
elseif fixed==1
    seldof=[seldof,4+[0.01,0.03]];
end
clear iSel4 

[T3,~,MasterDOF3]=tconstr(Constr,MasterDOF2);
clear Constr

% ## Apply boundary conditions
Constr=[zeros(size(seldof)).',ones(size(seldof)).',seldof.'];
[T4,~,MasterDOF4]=tconstr(Constr,MasterDOF3);
clear Constr

% ## SYSTEM REDUCTION
T=T1*T2*T3*T4;
Kr=T.'*K*T;
Mr=T.'*M*T;

%% ADD SPRINGS FOR THE ROLLER SUPPORTS AT SIDE LEUVEN

if fixed==0
   [~,I]=selectdof(MasterDOF4,[3.01;4.01]);
   for ind1=1:length(I)
      Kr(I(ind1),I(ind1))=Kr(I(ind1),I(ind1))+kxx; 
   end
end

end