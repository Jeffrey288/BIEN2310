function v0 = solarBVP_20851871(R0, tau)
% REMINDER: Rename your function to solarBVP_<8-digit-Student-ID>
% Do not change input arguments or return value 

GM = 2.96e-4; % (AU^3 / d^2)

% Use shooting method OR finite difference method. 
% Remember to set return value v0
% ----------------------------------
% Add code here
 
v0_guess = sqrt(GM/R0);
v0 = fsolve(@f, v0_guess);

    function diff = f(v)
        tau_res = solarIVP_20851871(R0, v, false);
        diff = tau_res - tau;
    end


% ----------------------------------

% Plot the orbit (you may call solveIVP to do so)
% ----------------------------------
% Add code here
 
 
solarIVP_20851871(R0, v0, true);
 


% ----------------------------------


end