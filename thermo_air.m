%{
Thermodynamic Calculations for Air at 
input Temperatures and Pressures
%}

function [Cp_air, rho_air, mu_air, k_air, beta_air] = thermo_air(T,P)
    %heat capacity of air at sea-level pressure and atmospheric temperature
    Cp_air = py.CoolProp.CoolProp.PropsSI('CP0MASS','P', P, 'T', T ,'air'); %J/(kg K)
    %density of air
    rho_air = py.CoolProp.CoolProp.PropsSI('D','P', P ,'T', T,'air'); %kg/m^3
    %air viscosity
    mu_air = py.CoolProp.CoolProp.PropsSI('V','P', P ,'T', T ,'air'); %Pa s | kg/(m s)
    %air thermal conductivity
    k_air =  py.CoolProp.CoolProp.PropsSI('CONDUCTIVITY' , 'P', P, 'T', T, 'air'); %W/(m K)
    %Thermal Expansion Coefficient
    beta_air = py.CoolProp.CoolProp.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT' , 'P', P , 'T', T, 'air'); %1/K
end
