% Example of a quasi-static trainpassage simulation using the finite 
% element model of railway bridge KW51 
%
% We consider the trainpassages and the model with parameters obtained 
% from the model updating described in the following reference, situation
% AFTER retrofit:
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

clc, clearvars, close all
% ## Add MATLAB functions for plotting
addpath('toolbox')

%% LOAD MODEL 

% ## Input parameters
Es_rel=1;       % Corresponding to a steel Young's modulus of 200 GPa
hb_rel=1.0833;  % Corresponding to a ballast thickness of 0.65 m     
kxx=0;          % Fixed boundary conditions, no springs are applied         
Id_rel=0.08;    % EIyy = 6.14*1e6 Nm^2 (out-of-plane bending stiffness)
Ad_rel=0.88;    % EA = 36.8*1e8 N (axial stiffness)   
fixed1=0;       % Support conditions as indicated on blueprints
fixed2=1;       % Hinges for all four supports
plot_elem=1;    % Plot nodes and elements
addmass=1;      % The boxes applied during the retrofit are modeled
diags_box=[];   % No individual boxes are added

% ## Load model with support conditions as indicated on blueprints
[Kr1,Mr1,T1,Nodes,Types,Sections,Materials,Elements,DOF]=...
    model_kw51(Es_rel,hb_rel,kxx,Id_rel,Ad_rel,fixed1,plot_elem,addmass,...
    diags_box);

% ## Load model with pinned supports
[Kr2,Mr2,T2,~,~,~,~,~,~]=...
    model_kw51(Es_rel,hb_rel,kxx,Id_rel,Ad_rel,fixed2,0,addmass,diags_box);

clear 'Es_rel' 'hb_rel' 'kxx' 'Id_rel' 'Ad_rel' 'fixed1' 'fixed2' 
clear 'plot_elem' 'addmass' 'diags_box'

%% LOAD SENSOR SELECTION MATRICES

% If not existing, compose sensor selection matrices
if exist('sensor_definition.mat','file')==0
    select_sensors
end

% Load sensor selection matrices
load('sensor_definition.mat')

%% TRAINPASSAGE SIMULATION

% ## Trainpassage info
% The following two passages are considered: 
% 1) Desiro MS08 train, track B, on 28 June 2020 at 06:50 UTC
% 2) Desiro MS08 train, track A, on 5 September 2020 at 21:11 UTC

% Timestamp assigned to the trainpassage
sdn_passage=[datenum(2020,06,28,06,50,26);datenum(2020,09,05,21,11,48)];
% Axle loads Desiro MS08 train in Newton ([1], table 3)
AL=[132.0,131.4,116.8,116.4,94.7,94.6,94.6,94.7,116.4,116.8,131.4,132.0]*1e3;
% Track: first track B ([1], figs. 8-10), next track A ([1], figs. 11-12)
tracks={'B';'A'};
% Direction: Leuven or Brussels
directions={'Leuven';'Brussels'};
% Axle distances with respect to axle 1
dxa=[0.0,2.3,18.6,20.9,26.1,28.4,44.7,47.0,52.2,54.5,70.8,73.1].';
% Axle location at the moment where the train enters the bridge (t=0 s)
% Different for both trainpassages
x0=[116.316+dxa,-1*dxa];
% Train speed at the moment when the train enters the bridge 
% (different for both passages, + = direction Leuven)
v0=[-27.95,27.84];
% Mean train acceleration 
% (different for both passages)
a0=[-0.343,-0.626];
% Construct time vector
t=-2:0.1:10;

% ## Simulate and plot the response
for indt=1:2 % Loop over the two passages

    track=tracks{indt};
    
    % Reconstruct axle location time history, measured curvilinearly along
    % the rails (xa = x0 + v0*t + a0*t^2/2, va = v0 + a0*t, aa = a0)
    xa=(x0(:,indt)+v0(indt)*t+a0(indt)*t.^2/2).';
    % Reconstruct train speed
    va=(v0(indt)+a0(indt)*t).';

    % Detect when first axle enters/all axles have left the bridge
    [~,i0]=min(abs(max(xa.')));       
    [~,i1]=min(abs(min(xa.')-116.316)); 

    % Visualize the axle locations
    newfig([20,10]) % Open new figure
    plot(xa,t,'b-','linewidth',0.5)
    xline([0,116.316],'k-')
    yline([t(i0),t(i1)],'k:')
    xlabel('Position [m]','fontsize',10)
    ylabel('Time [s]','fontsize',10) 
    xl=[-50,165]; xlim(xl);
    yl=get(gca,'YLim');
    text(-46,yl(2)-0.08*(yl(2)-yl(1)),...
         sprintf(' %s, track %s, direction %s, \n v0 = %0.2f km/h, a0 = %0.2f m/s^2, %i axles',...
           datestr(sdn_passage(indt),'dd/mm/yyyy HH:MM:SS'),...
           tracks{indt},directions{indt},v0(indt)*3.6,a0(indt),size(x0,1)));
    text(-3.5,yl(1)+0.015*(yl(2)-yl(1)),'Brussels','rotation',90,'fontsize',10)
    text(117.5,yl(1)+0.015*(yl(2)-yl(1)),'Leuven','rotation',90,'fontsize',10)

    % Determine the load vector for the obtained axle locations and speed
    P=xrail_2_P(track,xa,AL,va);
    
    % Solve the equilibrium equations

    % CASE 1: support conditions as indicated on blueprints
    % Reduce the load vector
    Pr1=T1.'*P;
    % Compute the quasi-static response (this may take a while...)
    Ur1=Kr1\Pr1;
    % Expand the displacement vector 
    U1=T1*Ur1;

    % CASE 2: pinned supports
    % Reduce the load vector
    Pr2=T2.'*P;
    % Compute the quasi-static response (this may take a while...)
    Ur2=Kr2\Pr2;
    % Expand the displacement vector 
    U2=T2*Ur2;

    % Plot the vertical deflection at mid-span for support cases 1 and 2
    x_sel=58.158; y_sel=0; z_sel=0.2719;
    % Select the corresponding node on girder A
    [~,iSel]=min(sqrt((Nodes(:,2)-x_sel).^2+(Nodes(:,3)-y_sel).^2+...
        (Nodes(:,4)-z_sel).^2));
    seldof=Nodes(iSel,1)+0.03; % Selected DOF (0.03 = vertical)
    clear iSel

    newfig % Open new figure
    plot(t,selectdof(DOF,seldof)*U1,'b--')
    hold on
    plot(t,selectdof(DOF,seldof)*U2,'r-.')
    xline([t(i0),t(i1)],'k:')
    xlim([-2,10]), set(gca,'xtick',-2:2:10)
    ylim([-6,1]*1e-3), set(gca,'ytick',(-6:1)*1e-3)
    xlabel('Time [s]')
    ylabel('Displacement [m]')
    % Add annotation
    annotation('TextBox','String','Mid-span displacement girder A', ...
        'Position',get(gca,'Position'),'LineStyle','none', ...
        'HorizontalAlignment','right');
    legend('Blueprint conditions','Pinned supports','location','southeast')

    % Compute and plot the corresponding strains, see [1] figs. 8 / 11
    chSel=1:8;
    strain1=Ss(chSel,:)*U1;
    strain2=Ss(chSel,:)*U2;
    labels_strain=labels_s(chSel);
    newfig([25,6]) % Open new figure
    for ind1=1:8
        subplot(2,4,ind1)
        plot(t,strain1(ind1,:),'b--')
        hold on
        plot(t,strain2(ind1,:),'r-.')
        yline(0,':')
        xline([t(i0),t(i1)],'k:')
        xlim([-2,10]), set(gca,'xtick',-2:2:10)
        ylim([-1,5]*1e-5), set(gca,'ytick',(-1:5)*1e-5)
        xlabel('Time [s]')
        ylabel('Strain [-]')
        % Add annotation
        annotation('TextBox','String',labels_strain{ind1}, ...
            'Position',get(gca,'Position'),'LineStyle','none', ...
            'HorizontalAlignment','right');
    end

    % Compute and plot the corresponding displacements at the supports
    % see [1] figs. 9 / 12
    chSel=1:2;
    displ1=Sd(chSel,:)*U1;
    displ2=Sd(chSel,:)*U2;
    labels_displ=labels_d(chSel);
    newfig([25,4.5]) % Open new figure
    for ind1=1:2
        subplot(1,2,ind1)
        plot(t,displ1(ind1,:)*1e3,'b--')
        hold on
        plot(t,displ2(ind1,:)*1e3,'r-.')
        yline(0,':')
        xline([t(i0),t(i1)],'k:')
        xlim([-2,10]), set(gca,'xtick',-2:2:10)
        ylim([-1,3]), set(gca,'ytick',(-1:3))
        xlabel('Time [s]')
        ylabel('Displacement [mm]')
        % Add annotation
        annotation('TextBox','String',labels_displ{ind1}, ...
            'Position',get(gca,'Position'),'LineStyle','none', ...
            'HorizontalAlignment','right');
    end

end