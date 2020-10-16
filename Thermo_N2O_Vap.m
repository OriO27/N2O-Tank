%{
Calculates and returns thermodynamic state information at an input
temperature and pressure. 
%}

function [Cp_vap,rho_vap,mu_vap,k_vap,beta_vap,enthalpy_vap,enthalpy_vap_sat,P_sat,MM,Z_vap,R_vap]=Thermo_N2O_Vap(T,P)

    %Heat Capacity
    Cp_vap = py.CoolProp.CoolProp.PropsSI('CP0MASS','T',T,'Q',1,'N2O'); %J/(kg K)
    %Density
    rho_vap = py.CoolProp.CoolProp.PropsSI('D','T', T,'Q', 1,'N2O');     %kg/m^3
    %Thermal Expansion Coefficient 
    beta_vap = py.CoolProp.CoolProp.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','T', T,'Q',1,'N2O'); %1/K 
    %Enthalpy 
    enthalpy_vap = py.CoolProp.CoolProp.PropsSI('H','T',T,'P',P,'N2O'); %J/kg 
    %Saturation Enthalpy 
    enthalpy_vap_sat = py.CoolProp.CoolProp.PropsSI('H','P',P,'Q',1,'N2O'); %J/kg 
    %Saturation Pressure 
    P_sat = py.CoolProp.CoolProp.PropsSI('P','T',T,'Q',1,'N2O'); %Pa 
    %Molar Mass 
    MM = py.CoolProp.CoolProp.PropsSI('M','T',T,'P',P,'N2O'); %kg/mol 
    %Compressibility Factor 
    Z_vap = py.CoolProp.CoolProp.PropsSI('Z','T',T,'P',P,'N2O'); 
    %Gas Constant 
    R_vap = py.CoolProp.CoolProp.PropsSI('GAS_CONSTANT','T',T,'P',P,'N2O'); %J/(mol K)  

    %%
    %The Following from... "Thermophysical properties of nitrous oxide" from IHS ESDU 

    %Viscosity from EQ 4.10 
    mu_vap = exp(3.3281+(-1.18237) * (1/(T/309.59) - 1)^(1/3)+(-0.055155)*(1/(T/309.59)-1)^(4/3)); %microNs/m^2 
    mu_vap = mu_vap/1000000; %Ns/m^2 

    %thermal conductivity 
    k_vap = exp(-7.0887+(-0.276962)*(1-(T/309.59))^(-2/3)+2.8872*(1-(T/309.59))^(-1/3)+16.6116* ... 
    (1-(T/309.59))^(1/3)+(-11.8221)*(1-(T/309.59))^(2/3)); %mW/m*K 
    k_vap = k_vap/1000; %W/(m K) 

end 