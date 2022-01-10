function [best_tau_fan, best_error] = solveTauFan(alpha, acceptable_error, tau_lambda, tau_ramp, tau_c, eta_mechanical, f, ...
    gama_nozzle, gama_turbine, gama_compressor, e_turbine, e_fan, PI_1, PI_2, A, B)
    % Holds the best Values achieve to be returned
    best_tau_fan = 0;
    best_error = 100;
    % initial values
    error_ = 100;
    loopIndex = 0;
    % initial Guess
    tau_fan = ((tau_lambda/tau_ramp) - tau_c + 1 + alpha)/(1+alpha);
    while (error_ > acceptable_error && loopIndex < 200)
        if(loopIndex > 0)
            tau_fan = tau_fan - 0.01;
            if(tau_fan <= 0)
                break
            end
        end
        % Compressor, Turbine & Fan Power Balance yields
        tau_turbine = 1 - ( (tau_ramp/tau_lambda)/ (eta_mechanical*(1+f)) ) ...
            *((tau_c-1) + (alpha*tau_fan-1));
        
        % subsitute PHI_1 of tau Turbine
        gama_nozzel_number = (gama_nozzle -1) / gama_nozzle;
        tau_turb_power     = (gama_turbine*(gama_nozzle-1)) / ...
            (gama_nozzle*(gama_turbine-1)*e_turbine);
        PHI_1 = ((-alpha *tau_ramp*A^0.5) / (2*tau_lambda*eta_mechanical*(1+f))) ...
            *( (PI_1^gama_nozzel_number * tau_turbine) - tau_turbine^(1-tau_turb_power) )^-1 ...
            *(PI_1^gama_nozzel_number - (1- tau_turb_power) * (tau_turbine^(-tau_turb_power)))^0.5;
        % subsitute PHI_2 of tau Fan
        gama_comp_number = (gama_compressor -1) / gama_compressor;
        PHI_2 = 0.5*B^0.5 * (PI_2^gama_comp_number - (1-e_fan) * tau_fan^(-e_fan)) ...
            / ( (PI_2^gama_comp_number*tau_fan) - tau_fan^(1-e_fan) )^0.5;
        
        % Condition of optimum fuel consumption 
        a1 = PHI_1 * (1+f) / alpha;
        a2 = -PHI_2;
        % Calculate Abs Error
        error_ = ( 1 - abs( min(a1, a2) / max(a1, a2) ) )* 100;
        if(error_ < best_error)
            % Update best_error to be returned
            best_error = error_;
            best_tau_fan = tau_fan;
        end
        
        % Increament Loop Indexer
        loopIndex = loopIndex + 1;
    end
end