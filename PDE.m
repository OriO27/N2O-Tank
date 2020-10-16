%{%Used to approximate the partial derivative [del(u)/del(rho)]_T 



function [del_u_rho] = PDE(T, rho, m) %Function is not running
    %%    
    % +- 1% of the density 
    if T >= 306 %sets upper limit for temperature
       T = 306; 
    end
    %%
    if rho <= 0 %return 0 #no approximate solution if density = 0 
        del_u_rho = 0;
    else
        %%
        % +- 1% of the density
        rho_plus_1 = 1.01 * rho; %kg/m^3 
        rho_minus_1 = 0.99 * rho; %kg/m^3 

        %using the high and low density values to calculate high and low  
        %internal energy values 
        U_high = py.CoolProp.CoolProp.PropsSI('U','T',T,'D',rho_plus_1,'N2O'); %J/kg 
        U_low = py.CoolProp.CoolProp.PropsSI('U','T',T,'D',rho_minus_1,'N2O'); %J/kg 

        del_u = (U_high-U_low); %J/kg | returns the specific change in internal energy  
        del_rho = (rho_plus_1-rho_minus_1); %kg/m^3 
        %%
        % approximate solution to the PDE 
        del_u_rho = (del_u/del_rho); %(J m^3)/(kg^2) 
    end
     
end

 