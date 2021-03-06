function output = polytropicEff(compressor_turbine, tau_pi, value, eta, gamma)
    
    if(compressor_turbine == "compressor")
        if(tau_pi == "tau")
            % Pi
            output = (eta*(value-1) + 1)^(gamma/(gamma-1));
        else
            % Tau
            output = (value^((gamma-1)/gamma) - 1)/eta + 1;
        end
    elseif(compressor_turbine == "turbine")
        if(tau_pi == "tau")
            % Pi
            output = ((value-1)/eta + 1)^(gamma/ (gamma-1));
        else
            % Tau
            output = eta*(value^((gamma-1)/gamma) - 1) + 1;
        end        
    end
end