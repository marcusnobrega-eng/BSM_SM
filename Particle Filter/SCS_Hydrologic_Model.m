% SCS-CN Non-Linear Reservoir MODEL
% Developer: Marcus Nobrega
%
% Description: This function solves the SCS-CN infiltration model coupled
% with the non-linear reservoir model approach used in the software SWMM.
% It assumes that the watershed is a plane with known width and length and
% rainfall is considered uniform in space in this plane. The excess of
% infiltration is propagated considering the watershed as a reservoir with
% flows only leaving the resevoir for depths larger than the initial
% abstraction h0_w. This reservoir, however, has a slope, which allows the
% calculation of a velocity using Manning's equation kinematic wave
% approximation.
%
% Rainfall File: A file called Rainfall_Intensity must be at the current matlab folder and
% must containt time(min) in the first column, and rainfall intensity
% (mm/h) in the second column. This .csv file will be read and converted
% into rainfall in the model. Moreover, rainfall must be equally spaced in
% time, such that a constant rainfall time-step must be input.
%
% Inputs:
%
% dt - Time-step [sec]
% routing_time - Simulation duration [min]
% width_w - Average watershed width [m]
% length_w - Average watershed length [m] A/width_w or vice-versa
% Lamba - SCS-CN Infiltration factor of Ia = Lambda * S
% CN_per - Average CN of pervious areas
% h0_w - Initial abstraction in [m]
% n_w - Manning's roughness coefficient in SI units
% Aimp - Rate of impervious connected areas to drainage area [0,1]
% slope_w - Outlet Slope in [m/m]
% baseflow - Constant baseflow [m3/s]
% ks - Saturated hydraulic conductivity [mm/hr]
% kr - Recovery rate or exfiltration rate in [mm/hr]
% A_GI - Ratio of GI areas (Area of GI / Total Watershed Area)
% CN_GI - Average CN of GI systems
% catch_GI_imp - Percentage of impervious areas that drains to GIs
%
% Outputs
% Graphs of flows and rainfall in .PDF
% Vectors of time [sec], flow discharge [m3/s], water depth [m], and
% infiltration rate [mm/h]
%
% Example:
% Calculate the hydrograph for a constant rainfall rate of 50 mm/h during
% 60 min. Assume ETP = 0.25 mm/h. Solve the SCS problem for a duration of
% 120 min. Use a time-step of 1 sec and the watershed properties are given
% as follows:
%
% width_w = [500]'
% length_w = [1000]'
% CN_per = [90]'
% h0_w = [0.01]'
% n_w = [0.02]'
% Aimp = [0.0]'
% slope_w = [0.015]'
% baseflow_w = [0.00]'; % (m3/s)
% ks_w = [50]';
% kr_w = [5]';
% A_GI = [0]'
% CN_GI = [65]'
% catch_GI_imp = [0]'
% flag_plot: 1 if you want to plot results
% Which results in:
%
% [time,Q_sum_watershed,H_w] =  SCS_Hydrologic_Model(1,120,500,500,90,0.01,0.02,0,0.015,0,50,5,0,65,0);

function [time,Q_watershed,H_w,f] = SCS_Hydrologic_Model(dt,routing_time,rainfall_data,time_step_rainfall,time_rainfall,ETP_table,time_step_ETP,time_ETP,width_w,length_w,Lambda,CN_per,h0_w,n_w,Aimp,slope_w,baseflow_w,ks_w,kr_w,A_GI,CN_GI,catch_GI_imp,flag_plot)
nw = 1; % Forcing the model only one time per watershed
%% Number of time-steps
nsteps = routing_time*60/dt;
%% Watershed Area
Area_w = width_w.*length_w; % watershed area in m2
%% Area Check
% Check catchment areas
if max(A_GI + Aimp) > 1
    error('Drainage area of impervious and GIs larger than the watershed')
end

%% Disagragete rainfall into model time-step
nsteps_rainfall = time_step_rainfall*(size(time_rainfall,1)-1)*60/dt; % I've changed it
rainfall_disagregated = zeros(nsteps_rainfall,nw);
% Disagregation
for i = 1:nsteps_rainfall
    z = ceil(i/(time_step_rainfall*60/dt));
    rainfall_disagregated(i,:) = rainfall_data(z,:);
end
rainfall_disagregated(nsteps_rainfall:nsteps,:) = repmat(zeros(1,nw),(nsteps - nsteps_rainfall + 1),1); % Increasing rainfall for time-steps larger than rainfall duration

%% Disagreagte ETP into model time-step
if isempty(ETP_table)
    ETP_disagregated = 0*rainfall_disagregated;
else
    % Disagregation
    nsteps_ETP = time_step_ETP*(size(time_ETP,1)-1)*60/dt;
    ETP_disagregated = zeros(nsteps_ETP,nw);
    for i = 1:nsteps_ETP
        z = ceil(i/(time_step_ETP*60/dt));
        ETP_disagregated(i,:) = ETP_table(z,:);
    end
    ETP_disagregated(nsteps_ETP:nsteps,:) = repmat(ETP_table(end,1),(nsteps - nsteps_ETP + 1),1); % Repeating ETP last value for time-steps larger than ETP duration
end
%% SCS Calculations - Watershed Hydrologic Model
% Empty duration (T = 4.5 / sqrt(Ks) )
empt_time = 4.5*1./(sqrt(ks_w/25.4)); % equation made in in/hr resulting in T in hours
% SCS Initial Abstraction Factor (i.e. 0.2*S = Lambda*S)
% Lambda = 0.2;
% GI Inflow Areas
A_per = 1 - A_GI - Aimp; % percentage of the watershed
f_IMP = Aimp;
f_GI = catch_GI_imp.*Aimp + A_GI; % Percentage of the watershed
f_inflow_GI = max(f_GI./A_GI,0); % Percentage of A_GI that drains to GI areas
% Pervious Areas
S_per = 25400./CN_per - 254; % Potential infiltration of pervious areas in (mm)
% GI Areas
S_GI = 25400./CN_GI - 254; % Potential infiltration of GI areas (mm)
% Accumulated Precipitation + Rainfall from Impervious Areas
Pac_GI = cumsum(f_inflow_GI'.*rainfall_disagregated*dt/3600); % Accumulated precipitation in (mm)
% Effective Cumulative Precipitation
Pef_GI = ((max(Pac_GI - Lambda*S_GI',0)).^2./(Pac_GI + (1-Lambda)*S_GI'))/1000; % Cumulative Effective precipitation in (m)
% Pef_GI(nsteps_rainfall:nsteps,:) = repmat(Pef_GI(nsteps_rainfall,:),(nsteps - nsteps_rainfall + 1),1); % Increasing rainfall for time-steps larger than rainfall duration
Pef_GI(nsteps_rainfall:nsteps,:) = 0; % Increasing rainfall for time-steps larger than rainfall duration

% Incremental Precipitation
dPef_GI(2:(size(Pef_GI,1)),1:size(Pef_GI,2)) = Pef_GI(2:end,:) - Pef_GI(1:(end-1),:); % Incremental precipitation for pervious areas per each time-step in (m)
dPef_GI(1,:) = 0; % Initial dPef is assumed as zero for the first time-step
%% Preallocating Arrays for Flows and Water Depths
% Flow Vectors
Q_watershed  = zeros(nsteps,nw);
Q = Q_watershed;
% Water Depths
H_w = zeros(nsteps,nw);
H_w_0 = zeros(nw,1);
T = zeros(1,nw); % Time with no rainfall in hr
f = zeros(nsteps,nw);
%% Watershed Routing
for i = 1:nsteps
    % Non-linear reservoir routing as SWMM model
    if i == 1
        P = zeros(nw,1); S = S_per; S1 = S; F = zeros(nw,1); % What about the initial conditions?
    else
        P = P + dt/3600*rainfall_disagregated(i,:)';
        F1 = P - P.^2./(P + S1); % cumulative infiltration in mm
        f(i,:) = (F1 - F)/(dt/3600); % infiltration rate in mm/hr
        % Refeshing states
        F = F1; S = max(0,S - f(i,:)'.*dt/3600); % mm

        % Checking Emptying time
        if min(rainfall_disagregated(i,:)) == 0 % At least one watershed has no rainfall
            for j = 1:nw % All watersheds
                z = double(rainfall_disagregated(i,nw) == 0); % watersheds with no rainfall
                T(j) = (T(j) + dt/3600)*z;
                f(i,j) = f(i-1,j);
                if T(j) > empt_time(j) % refresh values to original
                    P(j) = 0; F(j) = 0; S(j) = S(j); S1(j) = S(j); % Initial Values
                else
                    % First Order Recover
                    S(j) = S(j) + kr_w(j)*(S_per(j) - S(j))*dt/3600; % first order recovery rate
                end
            end
        end
    end
end
% Evaporation
evp = ETP_disagregated; % mm/hr
%% Watershed Routing
imp_flow = [rainfall_disagregated ; zeros(nsteps - size(rainfall_disagregated,1),nw)]/1000/3600'; % imp flow (m/s)
for i = 1:nsteps
    % Non-linear reservoir routing
    Q(i,1:nw) = 1./n_w'.*width_w'.*slope_w'.^(1/2).*(max(H_w_0' - h0_w',0)).^(5/3); % Flow in (m3/s)
    H_w(i,:) = H_w_0' + (rainfall_disagregated(i,:) - evp(i,:) - f(i,:))*dt/3600/1000 - Q(i,1:nw)./(Area_w'.*A_per')*dt + imp_flow(i,:)*dt.*f_IMP' + dPef_GI(i,:).*f_GI' ; % (m) % Mass balance with H_w in (m) considering imp areas, GI areas, and imp flows. Imp flow is routed with the non-linear reservoir
    % Refresh hw0
    H_w_0 = H_w(i,:)';
    if i < nsteps
        f(i+1,:) = min(f(i+1,:),rainfall_disagregated(i,:) + H_w_0'/dt*1000*3600); % mm/hr, available mass balance
    end
end
% Adding baseflow to the Watershed Flow.
% flow
Q(:,1:nw) = Q(:,1:nw) + baseflow_w'; %
% Defining Watershed Matrix
Q_watershed = Q(:,1:nw);
% Applying Mass balance at each node - Only Considering Flows from
% Watersheds
time = (dt:dt:routing_time*60)/60; % time in minutes

%% Plot Results
if flag_plot == 1
    clf
    shg
    colors = linspecer(3);
    set(gcf,'units','inches','position',[4,1,6.5,7])
    nplots = 3;

    % Plot Hydrographs for Watersheds
    if size(Q_watershed,2) > 0 % It means we have data
        subplot(nplots,1,1)
        plot(time,Q_watershed,'linewidth',3,'LineStyle','-','Color',colors(2,:))
        ylabel('Flow Discharge [m\textsuperscript{3}/s]','Interpreter','latex')
        hold on
        yyaxis("right")
        plot(time,rainfall_disagregated,'linewidth',2,'LineStyle','-.','Color',colors(1,:))
        ylabel('$i$ [mm/h]','Interpreter','latex')
        set(gca,'YDir','reverse');
        ylim([0 300]);
    end
    legend('Flow discharge','Rainfall Intensity','interpreter','latex','Location','best')
    xlabel('Elapsed time (min)','Interpreter','latex')
    set(gca,'ycolor','black')
    hold off

    if size(Q_watershed,2) > 0 % It means we have data
        subplot(nplots,1,2)
        plot(time,f,'linewidth',3,'LineStyle','-','Color',colors(2,:))
        ylabel('Infiltration Rate [mm/h]','Interpreter','latex')
        hold on
        plot(time,rainfall_disagregated,'linewidth',2,'LineStyle','-.','Color',colors(1,:))
        legend('Infiltration','Rainfall Intensity','interpreter','latex','Location','best')
        xlabel('Elapsed time (min)','Interpreter','latex')
        set(gca,'ycolor','black')
    end

    if size(Q_watershed,2) > 0 % It means we have data
        subplot(nplots,1,3)
        plot(time,(dt/3600)*cumsum(f),'linewidth',3,'LineStyle','-','Color',colors(3,:))
        ylabel('Volume [mm]','Interpreter','latex')
        hold on
        plot(time,(dt/3600)*cumsum(rainfall_disagregated),'linewidth',2,'LineStyle','-.','Color',colors(1,:))
        legend('Cumulative Infiltration','Cumulative Rainfall','interpreter','latex','Location','best')
        xlabel('Elapsed time (min)','Interpreter','latex')
        set(gca,'ycolor','black')
    end
    exportgraphics(gcf,'Results.pdf','ContentType','vector')
end
end