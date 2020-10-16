%{
SD_Main
Hybrid Oxidizer Filling Dynamics
N2O Tank Emptying
%}

clc; clear; clear all;
%%
%Tank Constraints
pyenv; %creates a python enviorment
global dVLiqhold; global LLiqhold

global g; g = 9.81; %m/s^2
%Tank Outter (ro) and Inner (ri) radii || Diameter/2 
global ro; ro = 0.1778/2; %m (7in)     
global ri; ri = 0.1524/2; %m (6in)

%Tank Properties (Material is Al for now)
global rho_wall; rho_wall = 2710; %kg/m3 | Density
global C_wall; C_wall = 887; %J/( kg K) | Heat Capacity
global k_wall; k_wall = 150; %W/(m K) | Thermal Conductivity
global L_tank; L_tank = 1.016; %m (40in) | Tank Length

%Calculates Volume of Tank Material 
global V_tank_wall; V_tank_wall = pi * (ro^2 - ri^2) * L_tank; %m^3

%Calculates Tank Holding Volume (Inner Volume)
global V_tank_hold; V_tank_hold = pi * ri^2 * L_tank; %m^3

%Atmospheric Conditions 
global T_atm; T_atm = 2 + 273.15; %K | Temperature
global P_atm; P_atm = 101325; %Pa | Pressure 
global P_out; P_out = 101325; %Pa | Pressure 
global P_tank; P_tank = 6894800; %Pa |Tank Pressure

%Initial Temperature Conditions
T_wall_vap = 20 + 273.15; %self.T_atm #K | wall temperature at vapor region
T_wall_liq = 20 + 273; %K | wall temperature at liquid region
T_lv_surf = -30 + 273.15; %K | inside temeprature at liquid vapor boundry
T_liq = 11 + 273.15; %K | liquid temperature
T_vap = 15 + 273.15; %K | vapor temperature

N2O_mass = 15.0; %kg | Total mass of N2O
fill_percent = 0.90; % | percentage initially filled with liquid N2O
t_total = 600; %s () | Total process time...

global M_dot;
%if Fill == True
M_dot =  N2O_mass / t_total;
%elseif Fill == False
%    M_dot = N2O_mass / t_total;
%end


%%
%%Initial Conditions
Tank_IC = [...
            0,... Q_liq_surf     
            0,... Q_surf_vap     
            0,... Q_wall_vap_out
            0,... Q_wall_liq_out     
            0,... m_evap
            0,... m_cond
            0,... m_vap
            0,... Q_in_vap
            0,... Q_in_liq
            py.CoolProp.CoolProp.PropsSI('D','Q',1,'P',P_out, 'N2O'),... rho_vap
            0,... V_wall_vap
            0,... Q_wall_vap_in
            0,... Q_wall_liq_in
            0,... Q_wall_vap_cond
            0,... Q_wall_liq_cond
            0,... m_wall_vap_in
            0,... m_wall_liq_in
            T_wall_vap,... T_wall_vap
            T_wall_liq,... T_wall_liq
            N2O_mass * fill_percent,... m_liq
            py.CoolProp.CoolProp.PropsSI('D','Q',0,'P',P_tank, 'N2O'),... rho_liq
            py.CoolProp.CoolProp.PropsSI('U','Q',1,'P',P_tank, 'N2O'),... U_vap
            py.CoolProp.CoolProp.PropsSI('U','Q',0,'P',P_tank, 'N2O'),... U_liq
            T_vap,... T_vap
            T_liq,... T_liq
            P_tank];% P_tank]

        
%%


global t_step; t_step = 0.1;

PlotInfo = Tank_IC;
dVLiqhold = [0, 0];
LLiqhold = [0,0];

time_span = [0];
time = 0;
global Counter
for Counter = 0:1000
    
    [dQ_liq_surf_dt, dQ_surf_vap_dt, dQ_wall_vap_out_dt, dQ_wall_liq_out_dt,...     
    dm_evap_dt, dm_cond_dt, dm_vap_dt, dQ_in_vap_dt, dQ_in_liq_dt, drho_vap_dt,...
    dV_wall_vap_dt, dQ_wall_vap_in_dt, dQ_wall_liq_in_dt, dQ_wall_vap_cond_dt,...
    dQ_wall_liq_cond_dt, dm_wall_vap_in_dt, dm_wall_liq_in_dt, dT_wall_vap_dt, dT_wall_liq_dt,...
    dm_liq_dt, drho_liq_dt, dU_vap_dt, dU_liq_dt, dT_vap_dt, dT_liq_dt, P_tank] = Tank_ODE_solve(Tank_IC);    
    
    Output = [dQ_liq_surf_dt, dQ_surf_vap_dt, dQ_wall_vap_out_dt, dQ_wall_liq_out_dt,...     
    dm_evap_dt, dm_cond_dt, dm_vap_dt, dQ_in_vap_dt, dQ_in_liq_dt, drho_vap_dt,...
    dV_wall_vap_dt, dQ_wall_vap_in_dt, dQ_wall_liq_in_dt, dQ_wall_vap_cond_dt,...
    dQ_wall_liq_cond_dt, dm_wall_vap_in_dt, dm_wall_liq_in_dt, dT_wall_vap_dt, dT_wall_liq_dt,...
    dm_liq_dt, drho_liq_dt, dU_vap_dt, dU_liq_dt, dT_vap_dt, dT_liq_dt, P_tank];
    

    NewOutput = Output.* t_step;
    NewOutput(1,26) = 0;
    
    Tank_IC = Tank_IC + NewOutput;
    disp(Tank_IC(1,21))
    disp(Tank_IC(1,25))
    TankPressure = py.CoolProp.CoolProp.PropsSI('P', 'D', Tank_IC(1,21), 'T', Tank_IC(1,25), 'N2O');
    if TankPressure < 90000
        TankPressure = 90000
    end
    disp(TankPressure)
    Tank_IC(1,26) = TankPressure
    PlotInfo = [PlotInfo ; Tank_IC];
    %disp('New Tank IC')
    %disp(Tank_IC)
    
    time = time + t_step;
    
    time_span = [time_span; time];

%%
%{
    figure
    plot(time_span, PlotInfo(:,1),'--g')
    xlabel('Time (s)')
    ylabel('QLiqSurf (J)')

    figure
    plot(time_span, PlotInfo(:,2),'r')
    xlabel('Time (s)')
    ylabel('QSurfVap (J)')
    
    figure
    plot(time_span, PlotInfo(:,3),'c')
    xlabel('Time (s)')
    ylabel('QWallVapOut (J)')
    %}
    figure
    plot(time_span, PlotInfo(:,24), 'b')
    hold on 
    plot(time_span, PlotInfo(:,25), 'r')
    title('Oxidizer Tank Temperature v Time')
    legend('T[Vap]','T[Liq]')
    xlabel('Time (s)')
    ylabel('Temperature (K)')
    
    
    figure
    plot(time_span, PlotInfo(:,7), 'c')
    hold on 
    plot(time_span, PlotInfo(:,20), 'm')
    title('Oxidizer Mass v Time')
    legend('m[Vap]','m[Liq]')
    xlabel('Time (s)')
    ylabel('Mass (kg)')
end