function output = adiabaticEff(compressor_turbine, tau_pi, value, e, gamma)
    
    if(compressor_turbine == "compressor")
        if(tau_pi == "tau")
            output = value^(gamma*e/(gamma-1));
        else
            output = value^((gamma-1)/(gamma*e));
        end
    elseif(compressor_turbine == "turbine")
        if(tau_pi == "tau")
            output = value^(gamma/(gamma-1)/e);
        else
            output = value^((gamma-1)*e/gamma);
        end        
    end
end