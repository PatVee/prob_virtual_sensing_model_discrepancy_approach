% Example of a modal analysis using the finite element model of railway 
% bridge KW51.
%
% We consider the model with parameters obtained from the model updating 
% described in the following reference, situation BEFORE retrofit:
%
% [1] K. Maes, L. Van Meerbeeck, E. Reynders, and G. Lombaert.  
%     Validation of vibration-based structural health monitoring on retro-
%     fitted railway bridge KW51. Mechanical Systems and Signal Processing, 
%     165(108380):1-24, 2022. https://doi.org/10.1016/j.ymssp.2021.108380.
%
% October 2025
%
% # AUTHORS
% Kristof Maes, Menno van de Velde and Geert Lombaert
% KU Leuven, Department of Civil Engineering, Structural Mechanics Section
%
% # UNITS 
% [kg, m, s]

close all, clear all, clc

% ## Add MATLAB functions for plotting
addpath('toolbox')

%% LOAD MODEL 

% ## Input parameters
Es_rel=1;       % Corresponding to a steel Young's modulus of 200 GPa
hb_rel=1.1;     % Corresponding to a ballast thickness of 0.66 m     
kxx=NaN;        % Fixed boundary conditions, no springs are applied         
Id_rel=0.06;    % EIyy = 4.6*1e6 Nm2 (out-of-plane bending stiffness)
Ad_rel=0.78;    % EA = 32.6*1e8 N (axial stiffness)   
fixed=1;        % Hinges for all four supports
plot_elem=1;    % Plot nodes and elements
addmass=0;      % The boxes applied during the retrofit are not modeled
diags_box=[];   % No individual boxes are added

% ## Load the model
[Kr,Mr,T,Nodes,Types,Sections,Materials,Elements,DOF]=...
    model_kw51(Es_rel,hb_rel,kxx,Id_rel,Ad_rel,fixed,plot_elem,addmass,...
    diags_box);
clear 'Es_rel' 'hb_rel' 'kxx' 'Id_rel' 'Ad_rel' 'fixed' 'plot_elem' 
clear 'addmass' 'diags_box'

%% MODAL ANALYSIS

% ## Compute the first 50 eigenmodes
[phir,omega]=eigfem(Kr,Mr,50);

% ## Expand the mode shapes to the original degrees of freedom
phi=T*phir;

% ## Plot the modes listed in table 6 of the MSSP reference 
mSel=[1:7,15,20,25,27,30,39,41];    % Mode numbers
fSel=omega(mSel)/2/pi;              % Natural frequencies [Hz]
newfig % Open new figure
for ind1=mSel(1:4) % Only the first 4 modes are plotted here
    subplot(2,2,ind1), plotdisp(Nodes,Elements,Types,DOF,phi(:,ind1),...
        'dispmax','off')
    title(sprintf('Mode %i - %0.2f Hz',ind1,omega(ind1)/2/pi))
end

% ## Extract the mode shapes at the sensor locations on the bridge deck
% and the arches (corresponding to the 'trackedmodes' published on Zenodo:
% https://doi.org/10.5281/zenodo.3745914).

% If not existing, compose sensor selection matrices
if exist('sensor_definition.mat','file')==0
    select_sensors
end

% Load sensor selection matrices
load('sensor_definition.mat')

% Select displacement mode shapes at the sensor locations
phia=Sa*phi(:,mSel);
