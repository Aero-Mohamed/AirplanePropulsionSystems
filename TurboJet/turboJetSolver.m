%% Clear Workspace
% Mohamed Ahmed Hassan Ahmed
clear
clc

%% Defintion of Aussumption variable
% value from state of the art table
% Level of Technology 4
pi_diff             = 0.96;
pi_nozzle           = 0.97;
pi_fan_diff         = 0.995;
e_compressor        = 0.9;
e_fan               = 0.89;
pi_burner           = 0.96;
eta_burner          = 0.999;
e_turbine           = 0.9;
eta_mechanical      = 0.995;

gama_nozzle         = 1.4;
gama_compressor     = 1.4;
gama_turbine        = 1.33;

R                   = 287;
h_PR                = 4.28e7;
nozzleConfiguration = "C"; % "C" or "CD"

%% Defintion of Givens variables
% Operating Condition
T_0                 = 288; %kelvin @ 32000 ft
P_0                 = 101325; % pa @ 32000 ft
a_0                 = 340.17; % m/s
M_0                 = 0.0;
V_0                 = M_0 * a_0;
% Design Parameter
pi_compressor       = 6.458;
T_t4                = 1170; % kelvin

%% Defintion of Variables to be calculated from givens & assumption

Cp_t                = gama_turbine/(gama_turbine-1) * R;
Cp_c                = gama_compressor/(gama_compressor-1) * R;
tau_ramp            = (1 + ( (gama_nozzle-1)/2*M_0^2) );
pi_ramp             = polytropicEff("compressor", "tau", tau_ramp, 1, gama_compressor); % Eff
tau_lambda          = (Cp_t*T_t4)/(Cp_c*T_0);
tau_c               = adiabaticEff("compressor", "pi", pi_compressor, 0.88, gama_compressor); % Eff
% First Law of thermodynamics across Burner yields
f                   = ((tau_c*tau_ramp) - tau_lambda) / (tau_lambda - (eta_burner*h_PR/(Cp_c*T_0)));
% Turbine Compressor Power Balance
tau_t               = 1 - tau_ramp*(tau_c-1)/(eta_mechanical*(1+f)*tau_lambda);
pi_turbine          = polytropicEff("turbine", "tau", tau_t, 0.9, gama_turbine); % Eff


%% Checking The Nozzle 
Pt9_P0 = pi_nozzle*pi_turbine*pi_burner*pi_compressor*pi_diff*pi_ramp;
if(nozzleConfiguration == "C")
    % Assuming Convergent Nozzle
    if( Pt9_P0 > 1.89 )
        disp("Chocked Nozzle Pt9_P0 = "+ char(vpa(Pt9_P0, 6)) +" > 1.89")
        M_9     = 1;
        Pt9_P9  = ( 1 + ((gama_turbine-1)/2))^(gama_turbine/(gama_turbine-1));
        P9_P0   = Pt9_P0 / Pt9_P9;
    elseif (Pt9_P0 == 1.89)
        disp('Chocked Nozzle Pt9_P0 = 1.89')
        M_9     = 1;
        P9      = P_0;
        P9_P0   = 1;
        Pt9_P9  = Pt9_P0;
    elseif( Pt9_P0 < 1.89 )
        disp('Unchocked Nozzle P9 = P0')
        P9      = P_0;
        P9_P0   = 1;
        Pt9_P9  = Pt9_P0;
        M_9     = sqrt(2*(Pt9_P9^((gama_turbine-1)/gama_turbine) - 1)/(gama_turbine-1));
    end
else
    % Assuming Convergent-Divargent Nozzle & Choked.
    disp('Choked C-D-Nozzle P9 = P0 and M_9 > 1')
    P9      = P_0;
    P9_P0   = 1;
    Pt9_P9  = Pt9_P0;
    M_9     = sqrt(2*(Pt9_P9^((gama_turbine-1)/gama_turbine) - 1)/(gama_turbine-1));
end
T_t5 = T_t4*tau_t;
T_9 = (T_t5) / (Pt9_P9)^((gama_turbine-1)/gama_turbine);
a_9 = sqrt(gama_turbine * R * T_9);
V9 = M_9 * a_9;


%% Outputing S.Thrust, SFC & efficiencies
F_sp            = a_0 * ((1+f) * (V9/a_0) - M_0 + (1+f)*(T_9/T_0)*(1-1/P9_P0)/gama_compressor/(V9/a_0) );
SFC             =  f / F_sp;
if(V_0 > 0)
    eta_total       = V_0 / SFC / h_PR;
    eta_propulsive  = 2 / ((F_sp/V_0) + 2);
    eta_thermal     = eta_total / eta_propulsive;
else
    eta_thermal     = (F_sp + 2*V_0) / (2*SFC*h_PR);
    eta_total       = 0;
    eta_propulsive  = 0;
end
disp("Specific Thrust = " + F_sp);
disp("Specific Fuel Consumption (SFC) = " + char(vpa(SFC, 6)));
disp("Total Efficiency (eta_0) = " + char(vpa(eta_total, 6)));
disp("Propulsive Efficiency (eta_p) = " + char(vpa(eta_propulsive, 6)));
disp("Thermal Efficiency (eta_th) = " + char(vpa(eta_thermal, 6)));

