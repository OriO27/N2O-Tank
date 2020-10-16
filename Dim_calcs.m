%{
Function to calculate Tank dimensions

%}
function [V_liq, V_vap, L_liq, L_vap, A_liq_in, A_vap_in, A_liq_out,...
          A_vap_out, m_wall_liq, m_wall_vap, V_wall_vap, V_wall_liq]...
          = Dim_calcs(m_liq, rho_liq)
    
    global V_tank_hold
    global L_tank
    global ri
    global ro
    global rho_wall
    
    V_liq = m_liq / rho_liq; %m^3 | volume of liquid in tank
    if V_liq >= V_tank_hold
        V_liq = V_tank_hold;
    end

    V_vap = V_tank_hold - V_liq; %m^3 | Constrained volume of vapor in tank
    
    if V_liq + V_vap > V_tank_hold
        disp("Error in volume calcs...")
    end
    
    L_liq = V_liq / (pi * ri^2); %m | Length of liquid region
    if L_liq >= L_tank
        L_liq = L_tank;
    end
    
    
    L_vap = L_tank -L_liq ; %m | Constrained length of vapor region

    if L_liq + L_vap > L_tank
        disp("Error in length calcs...")
    end

    %Note::: Surface area does not include end caps; hollow tube assumed
    A_liq_in = 2 * pi * ri * L_liq; %m^2 |inner surface area of tank liquid region
    A_vap_in = 2 * pi * ri * L_vap; %m^2 |inner surface area of tank vapor region

    A_liq_out = 2 * pi * ro * L_liq; %m^2 | outter surface area of tank liquid region
    A_vap_out = 2 * pi * ro * L_vap; %m^2 | outter surface area of tank vapor region

    m_wall_liq = pi * (ro^2 - ri^2) * L_liq * rho_wall; %kg | mass of wall in liquid region
    m_wall_vap = pi * (ro^2 - ri^2) * L_vap * rho_wall; %kg | mass of wall in liquid region

    V_wall_vap = pi * (ro^2 - ri^2) * L_vap; %m^3 | volume of vapor region wall
    V_wall_liq = pi * (ro^2 - ri^2) * L_liq; %m^3 | volume of liquid region wall


end