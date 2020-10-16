%{
Computes ODE
%}
%%
function [dQ_liq_surf_dt, dQ_surf_vap_dt, dQ_wall_vap_out_dt, dQ_wall_liq_out_dt,...     
    dm_evap_dt, dm_cond_dt, dm_vap_dt, dQ_in_vap_dt, dQ_in_liq_dt, drho_vap_dt,...
    dV_wall_vap_dt, dQ_wall_vap_in_dt, dQ_wall_liq_in_dt, dQ_wall_vap_cond_dt,...
    dQ_wall_liq_cond_dt, dm_wall_vap_in_dt, dm_wall_liq_in_dt, dT_wall_vap_dt, dT_wall_liq_dt,...
    dm_liq_dt, drho_liq_dt, dU_vap_dt, dU_liq_dt, dT_vap_dt, dT_liq_dt, P_tank] = ...
    Tank_ODE_solve(Input)

 InputCell = num2cell(Input);
 %%
    global ri; global ro; global T_atm; global k_wall; global rho_wall; global C_wall;
    global M_dot; global P_tank; global g; global P_atm; global t_step; global dVLiqhold; global LLiqhold
    global Counter; global dVVaphold;

[Q_liq_surf, Q_surf_vap, Q_wall_vap_out, Q_wall_liq_out,...     
    m_evap, m_cond, m_vap, Q_in_vap, Q_in_liq, rho_vap,V_wall_vap, Q_wall_vap_in, Q_wall_liq_in,...
    Q_wall_vap_cond, Q_wall_liq_cond, m_wall_vap_in, m_wall_liq_in, T_wall_vap, T_wall_liq,...
    m_liq, rho_liq, U_vap, U_liq, T_vap, T_liq, P_tank] = InputCell{:} ;




    A_boundry = pi * (ri)^2; % surface area of the liquid vapor boundry

    [V_liq, V_vap, L_liq, L_vap, A_liq_in, A_vap_in, A_liq_out,...
    A_vap_out, m_wall_liq, m_wall_vap, ~, ~] = ...
    Dim_calcs(m_liq, rho_liq);
       
    %m^3/s
    dVLiqhold = [dVLiqhold, V_liq];    
    if Counter == 0
        dV_liq = 0;
    else
        dV_liq = (dVLiqhold(1,end) - dVLiqhold(1,end-1))/t_step;
    end
    
    dVVaphold = [dVVaphold, V_liq];
    if Counter == 0
        dV_vap = 0;
    else
        dV_vap = (dVVaphold(1,end) - dVVaphold(1,end-1))/t_step;
    end


    %m/s
    LLiqhold = [dVLiqhold, L_liq];
    if Counter == 0
        dL_liq_dt = 0;
    else
        dL_liq_dt = (LLiqhold(1,end) - LLiqhold(1,end-1))/t_step;
    end
    
    dL_vap_dt = - dL_liq_dt;

    %%
    %Temperature Corrections...
    if T_liq >= 309.4
        T_liq = 309;
    elseif T_liq <= 182.5
        T_liq = 183;
    end
    if T_vap >= 309.4
        T_vap = 309;
    elseif T_vap <= 182.5
        T_vap = 185;
    end
    
    %%
    
    if T_wall_liq >= 309.4
        T_wall_liq = 309;
    elseif T_wall_liq <= 182.5
        T_wall_liq = 183;
    end
    if T_wall_vap >= 309.4
        T_wall_vap = 309;
    elseif T_wall_vap <= 182.5
        T_wall_vap = 185;
    end
    
    %%
    %film temperature of tank wall and atmosphere 
    T_film_wall_vap_in = 0.5 * (T_wall_vap + T_atm);
    
    T_film_wall_liq_in = 0.5 * (T_wall_liq + T_atm);
    
    
    if T_film_wall_vap_in <= 185 %hacky approach to fix initial temperature starting at 0K due to evacuated tank... 
        T_film_wall_vap_in = 185; %%% NEEDS TO BE CHANGED TO T_OUT (Temperature at exit...)
    end
    %%
    %film temperature of liq and vap regions
    T_film_liq_surf = 0.5 * (T_vap + T_liq);
    if T_film_liq_surf <= 183
        T_film_liq_surf = 185; 
    end

    %%
    %HT Coefficients

    %Liq HT Coeffs 
    [Cp_liq, ~, mu_liq, k_liq, beta_liq, enthalpy_liq, enthalpy_liq_sat] = Thermo_N2O_Liq(T_film_liq_surf, P_tank);
    %Vap HT Coeffs
    [Cp_vap, rho_vap,mu_vap, k_vap, beta_vap, ~, enthalpy_vap_sat, ~, ~, ~, ~] = Thermo_N2O_Vap(T_film_liq_surf, P_tank);
    
    N2O_latent_heat = enthalpy_vap_sat - enthalpy_liq_sat;

    %%
    %Begin EQ Solving
    %%
    %%
    %1 J/s
    global E 
    if L_liq <= 0
        dQ_liq_surf_dt = 0;
    else
        dQ_liq_surf_dt = (0.15 * ((Cp_liq * rho_liq^2 * g * beta_liq * abs(T_film_liq_surf - T_liq ) * L_liq^3 )...
                / (mu_liq * k_liq ) )^(1/3) * (k_liq/L_liq)) * A_boundry * (T_film_liq_surf - T_liq);
    end
    
    %%
    %2 J/s
    if L_vap <= 0
        dQ_surf_vap_dt = 0;
    else
        dQ_surf_vap_dt = (0.15 * ((Cp_vap * rho_vap^2 * g * beta_vap * abs(T_film_liq_surf - T_vap) * L_vap^3 )...
            / (mu_vap * k_vap ) )^(1/3) * (k_vap/L_vap)) * A_boundry * (T_film_liq_surf - T_vap);
    end
    
    %%
    %3 J/s
    T_film_wall_vap_out = 0.5 * (T_vap + T_wall_vap);
    if T_film_wall_vap_out >= 309.4
        T_film_wall_vap_out = 309;
    elseif T_film_wall_vap_out <= 182.5
        T_film_wall_vap_out = 183;
    end
        
    [Cp_vap_w,~,mu_vap_w,k_vap_w,beta_vap_w,~,~,~,~,~,~]= Thermo_N2O_Vap(T_film_wall_vap_out, P_tank);
    if L_vap <= 0
        dQ_wall_vap_out_dt = 0;        
    else
        dQ_wall_vap_out_dt = (0.021 * ((Cp_vap_w * rho_vap^2 * g * beta_vap_w * abs(T_film_wall_vap_out - T_vap) * L_vap^3  )...
                / (mu_vap_w * k_vap_w ) )^(2/5) * (k_vap_w/L_vap)) * A_vap_in * (T_film_wall_vap_out - T_vap);
    end

    %%
    %4 J/s
    T_film_wall_liq_out = 0.5 * (T_liq + T_wall_liq);
    if T_film_wall_liq_out >= 309.4
        T_film_wall_liq_out = 309;
    elseif T_film_wall_vap_out <= 182.5
        T_film_wall_liq_out = 183;
    end
    [Cp_liq_w, ~, mu_liq_w, k_liq_w, beta_liq_w, ~, ~] = Thermo_N2O_Liq(T_film_wall_liq_out, P_tank);
    if L_liq <= 0
        dQ_wall_liq_out_dt = 0;
    else
        dQ_wall_liq_out_dt = (0.021 * ((Cp_liq_w * rho_liq^2 * g * beta_liq_w * abs(T_film_wall_liq_out - T_liq) * L_liq^3  )...
            / (mu_liq_w * k_liq_w ) )^(2/5) * (k_liq_w/L_liq)) * A_liq_in * (T_film_wall_liq_out - T_liq);
    end
    
    %%
    %5 kg/s
   
    dm_evap_dt = (dQ_liq_surf_dt - dQ_surf_vap_dt) / (N2O_latent_heat + (enthalpy_liq_sat - enthalpy_liq) );
    
    %%
    %6 kg/s
    
    [~,~,~,~,~,~,~,P_sat,MM,Z_vap,R_vap] = Thermo_N2O_Vap(T_vap, P_tank);
    if m_vap <= 0
        dm_cond_dt = 0;
    elseif P_tank > P_sat
        dm_cond_dt = ( ( P_tank - P_sat) * V_vap * MM ) / (Z_vap * R_vap  * T_vap * t_step );
    elseif P_tank <= P_sat 
        dm_cond_dt = 0;
    end

   %%
   %7 kg/s
   dm_vap_dt = dm_evap_dt - dm_cond_dt;
   %returns 0 if there is no vapor to be condensed AND the calculation
   %calls for a negative rate of change (i.e. turning into liquid)
   if m_vap < 0 & dm_vap_dt < 0
      dm_vap_dt = 0; 
   end
   
   %%
   %8 J/s
   dQ_in_vap_dt = dQ_wall_vap_out_dt + dQ_surf_vap_dt ;
   
   %%
   %9 J/s
   dQ_in_liq_dt = dQ_wall_liq_out_dt + dQ_liq_surf_dt;
   
   %%
   %10 kg/(m^3 s)
   if V_vap <= 0
       drho_vap_dt = 0;
   elseif rho_vap <= 0
        drho_vap_dt = 0;
   else
        drho_vap_dt = 1/V_vap * dm_vap_dt - m_vap / (V_vap)^2 * dV_vap;
   end
   

   %%
   %11 m^3/s
   dV_wall_vap_dt = pi * (ro^2 - ri^2) * dL_vap_dt;
   
   %%
   %12 J/s

   [Cp_air, rho_air, mu_air, k_air, beta_air] = thermo_air(T_film_wall_vap_in, P_atm);
   
   
   
   if L_vap <= 0
    dQ_wall_vap_in_dt = 0;
   else
    dQ_wall_vap_in_dt = (0.59 * (  (Cp_air * (rho_air )^2 * g * beta_air * abs(T_wall_vap - T_atm) * L_vap^3)...
        / (mu_air * k_air ) )^(0.25) * (k_air/L_vap)) * A_vap_out * (T_wall_vap - T_atm);
   end

    %%
    %13 J/s
  [Cp_air_l, rho_air_l, mu_air_l, k_air_l, beta_air_l] = thermo_air(T_film_wall_liq_in, P_atm);
  if L_liq <= 0
        dQ_wall_liq_in_dt = 0;
  else
        dQ_wall_liq_in_dt = (0.59 * (  (Cp_air_l * (rho_air_l )^2 * g * beta_air_l * abs(T_wall_liq - T_atm) * L_liq^3)... 
            / (mu_air_l * k_air_l ) )^(0.25) * (k_air_l/L_liq)) * A_liq_out * (T_wall_liq - T_atm);
  end

  %%
  %14 J/s

    if L_vap <= 0 && L_liq == 0
        dQ_wall_vap_cond_dt = 0;
    else
        dQ_wall_vap_cond_dt = (k_wall * (T_wall_liq - T_wall_vap) * pi * (ro^2 - ri^2))... 
            / (0.5 * L_liq + 0.5 * L_vap);
    end
    
    %%
    %15 J/s
    if L_vap <= 0 && L_liq == 0
        dQ_wall_liq_cond_dt = 0;
    else
        dQ_wall_liq_cond_dt = (k_wall * (T_wall_vap - T_wall_liq) * pi * (ro^2 - ri^2))...
            / (0.5 * L_liq + 0.5 * L_vap);
    end
   
    %%
    %16 kg/s
    dm_wall_vap_in_dt = dV_wall_vap_dt * rho_wall;
    
    %%
    %17 kg/s
    dm_wall_liq_in_dt = - dm_wall_vap_in_dt;
    
    %%
    %18 K/s
    if m_wall_vap <= 0
        dT_wall_vap_dt = T_wall_vap;
    else

        dT_wall_vap_dt = (dQ_wall_vap_in_dt - dQ_wall_vap_out_dt + dQ_wall_vap_cond_dt...
            + dm_wall_vap_in_dt * C_wall * (T_wall_liq - T_wall_vap))...
            / (m_wall_vap * C_wall);
    end
    
    %%
    %19 K/s
    if m_wall_liq <= 0
        dT_wall_liq_dt = T_wall_liq;
    else
        dT_wall_liq_dt = (dQ_wall_liq_in_dt - dQ_wall_liq_out_dt + dQ_wall_liq_cond_dt...
            + dm_wall_liq_in_dt * C_wall * (T_wall_vap - T_wall_liq))...
            / (m_wall_liq * C_wall);
    end
    
    %%
    %20 kg/s
    dm_liq_dt = - dm_evap_dt + dm_cond_dt - M_dot;
    
    %%
    %21 kg/(m^3 s)

    if V_liq <= 0 %|temporarily fixes no liquid volume .....
        drho_liq_dt  = rho_liq;
    else
        
        drho_liq_dt  = 1/V_liq * dm_liq_dt - (m_liq / V_liq^2) * (dV_liq);
        P_correct = 0.9994 * P_tank; %corrects coolprop +/- 1e-4% saturation proximity error
        PropCheck = py.CoolProp.CoolProp.PropsSI('D', 'P', P_correct, 'T', T_liq,'N2O');
        if (t_step * drho_liq_dt + rho_liq) >= PropCheck
            drho_liq_dt = 0;
        end
    end

    %%
    %22 J/s
    [h_evap,h_cond,h_out] = Enthalpy_calc(P_tank);
    dU_vap_dt  = dm_evap_dt * h_evap - dm_cond_dt * h_cond - P_tank * dV_vap + dQ_in_vap_dt;
    %%
    %23 J/s
    dU_liq_dt  = M_dot * h_out - dm_evap_dt * h_evap + dm_cond_dt * h_cond - P_tank * dV_liq + dQ_in_liq_dt;
    
    %%
    %24
    
    if m_vap <= 0 % return 0 for temp if theres no N2O...
        dT_vap_dt = 0;
    else
        dT_vap_dt = (1 / Cp_vap) * (1/m_vap * (dU_vap_dt - (U_vap/m_vap) * dm_vap_dt )...
            - drho_vap_dt * PDE(T_vap, rho_vap, m_vap) );
    end 
    
    %%
    %25
    
    if m_liq <= 0 % return 0 for temp if theres no N2O...
        dT_liq_dt  = 0;
    else
        dT_liq_dt = (1 / Cp_liq) * (1/m_liq * (dU_liq_dt - (U_liq/m_liq) * dm_liq_dt )...
            - drho_liq_dt * PDE(T_liq, rho_liq, m_liq) );
    end 

    
    %%
    %26
    P_tank = P_tank;
    

end