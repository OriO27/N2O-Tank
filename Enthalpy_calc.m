%{
calculates liquid and vapor enthalpy and input pressure
%}

function [h_evap,h_cond,h_out] = Enthalpy_calc(P)
    %%
    %Should eventually vary with Pressure instead of density.... 

    h_evap = py.CoolProp.CoolProp.PropsSI('H','P',P,'Q',1,'N2O'); %J/kg | evaporation specific enthalpy 
    h_cond = py.CoolProp.CoolProp.PropsSI('H','P',P,'Q',0,'N2O'); %J/kg | condensation specific enthalpy 

    %h_out ...will be intresting to model with inward flow since it will be dependent on the source  
    %tank exit conditions and the losses which occur in the line.  
    h_out = py.CoolProp.CoolProp.PropsSI('H','P',P,'Q',0,'N2O'); %J/kg  | outlet specific enthalpy | Assumes liquid is going in/out

end 