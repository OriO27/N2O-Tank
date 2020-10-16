%{
Calculates and returns thermodynamic state information of an input
temperature and pressure    
PROPSSI Solving as if quality (x=0) for liquid properties... may not be accurate
will also try having T as an input, or to calculate with P input and T input and average both values
%}

function [Cp_liq, rho_liq, mu_liq, k_liq, beta_liq, enthalpy_liq, enthalpy_liq_sat] = Thermo_N2O_Liq(T, P)

    %Heat Capacity
    Cp_liq = py.CoolProp.CoolProp.PropsSI('CP0MASS','T', T,'Q', 0 ,'N2O'); %J/(kg K)
    %density 
    rho_liq = py.CoolProp.CoolProp.PropsSI('D','T', T ,'Q', 0,'N2O'); %kg/m^3
    %Thermal Expansion Coefficient
    beta_liq = py.CoolProp.CoolProp.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT' , 'T', T , 'Q', 0, 'N2O'); %1/K
    %Enthalpy
    enthalpy_liq = py.CoolProp.CoolProp.PropsSI('H', 'T', T, 'P', P, 'N2O'); %J/kg
    %Saturation Enthalpy
    enthalpy_liq_sat = py.CoolProp.CoolProp.PropsSI('H', 'P', P, 'Q', 0, 'N2O'); %J/kg 
    
    %enthalpy_liq1 = PropsSI('H', 'T', T, 'P', P, 'N2O'); %J/kg
    %Saturation Enthalpy
    %enthalpy_liq_sat1 = PropsSI('H', 'P', P, 'Q', 0, 'N2O'); %J/kg 
    
    %P_sat1 = PropsSI('P', 'T', T, 'Q', 0, 'N2O'); %J/kg 

    %T_sat1 = PropsSI('T', 'P', P, 'Q', 0, 'N2O'); %J/kg 

    %%
    %The Following from... "Thermophysical properties of nitrous oxide" from IHS ESDU
    %Viscosity from EQ 4.9
    theta = ((309.59) -(5.24))/(T  -(5.24) );
    mu_liq = 0.0293423 * exp((1.6089) * (theta - 1)^(1/3) + 2.0439*(theta-1)^(4/3)); %mNs/m^2
    mu_liq = mu_liq/1000; %Ns/m^2
    % Thermal conductivity from EQ 4.11 
    % Note:: Literature suggests temp range ends at 10C, model below ends at T_crit... 
    % May Cause Instability in response.... 
    k_liq = 72.35 * (1 + 1.5 * (1 - (T/309.59))^(1/3) + (-3.5) * (1-(T/309.59))^(2/3) + 4.5 * (1- (T/309.59))); %W/(m K)
    k_liq = k_liq/1000; %W/(m K)

    
end