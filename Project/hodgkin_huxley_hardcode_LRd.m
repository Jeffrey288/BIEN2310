function res = hodgkin_huxley_hardcode_LRd(tf, I_t, params, options)
% Source: https://www.ahajournals.org/doi/10.1161/01.res.68.6.1501?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed
% http://www.scholarpedia.org/article/Models_of_cardiac_cell#fig:LRd_0d.gif
% Some reference is taken from
% https://github.com/opencor/speed-comparison/blob/master/luo_rudy_1991.m

    arguments
        tf
        I_t
        params
        options.Vm0 (1,1) double = NaN
        options.accuracy (1,1) int16 = 0
        options.getFreq (1,1) int8 = false
        % if getFreq is on, will only return frequency
    end

%   For this model, the time scale is ms
    initX = zeros(9,1);
    initX(1) = -85;
    initX(2) = 4e-4;
    initX = fsolve(@(x) simX(0, x), initX);
    initX(1) = -85;
    initX(2) = 4e-4;

    if (options.accuracy == 0)
        opt = odeset('MaxStep', 0.3);
    elseif (options.accuracy == 1)
        opt = odeset('RelTol', 1e-5, 'MaxStep', 0.2);
    else
        opt = odeset('RelTol', 1e-7, 'MaxStep', 0.1);
    end
    [tt, XX] = ode15s(@simX, [0, tf], initX, opt);

    % find the period
%     ind = islocalmax(XX(:, 1)) .* (XX(:,1) > 0); % element-wise and
    ind = islocalmax(XX(:, 1)) .* (XX(:,1) > 15); % element-wise and
    t_spikes = tt(logical(ind));
    if (length(t_spikes) <= 1)
        period = NaN;
    else
        l = length(t_spikes);
        period = mean(t_spikes(2:l) - t_spikes(1:(l-1)));
    end
    
    if (options.getFreq)
        res = 1 / period * 1000;
        return;
    end

%     figure(1)
%     plot(tt, arrayfun(I_t, tt), "c-");
%     title("Injected Current", "Interpreter", "latex")
%     xlabel("time (ms)")
%     ylabel("injected current (A/cm^2)")

    figure(1)
    plot(tt, XX(:, 1), "b-")
    title("Original HH: Membrane Voltage", "Interpreter", "latex")
    subtitle(sprintf("spike period = %.1f ms", period), "Interpreter", "latex")
    xlabel("time (ms)")
    ylabel("membrane potential difference (mV)")
    ylim([-30, 130])

    function XXdot = simX(t, XX)

        % gating variables
        Vm = XX(1);
        Ca_conc = XX(2);
        h = XX(3);
        jj = XX(4);
        m = XX(5);
        d = XX(6);
        f = XX(7);
        X = XX(8);
        Ki = XX(9);

        VNa = 54.4; % mV
        K_conc = 5.4; % mM
        VK = -77; % mV
        Cm = 1; %uF/cm^2
        VKi = (8.31*(37+273)/9.64e4)*log(5.4/145) * 1000; % mV
        % all currents are in uA/cm^2

        % Fast Sodium Current in Phase 1
        % m, h: fast activation and inactivation
        % jj: slow inactivation gate

        INa = 23 * m^3 * h * jj * (Vm - VNa);

        if (Vm > -40)
            alpha_h = 0;
            beta_h = 1/(0.13*(1+exp((Vm+10.66)/-11.1)));
            alpha_jj = 0;
            beta_jj = 0.3*exp(-2.535e-7*Vm)/(1+exp(-0.1*(Vm+32)));
        else
            alpha_h = 0.135*exp((80+Vm)/-6.8);
            beta_h = 3.56*exp(0.079*Vm)+3.1e5*exp(0.35*Vm);
            alpha_jj = (-1.2714e5*exp(0.2444*Vm)-3.474e-5*exp(-0.04391*Vm))*(Vm+37.78)/(1+exp(0.311*(Vm+79.23)));
            beta_jj = 0.1212*exp(-0.01052*Vm)/(1+exp(-0.1378*(Vm+40.14)));
        end
        alpha_m = 0.32*(Vm+47.13)/(1-exp(-0.1*(Vm+47.13)));
        beta_m = 0.08*exp(-Vm/11);

        % Slow inward current due to calclium inflow in phase 2 and 3
        % d, f: activation and inactivation
        
        Vsi = 7.7 - 13.0287 * log(Ca_conc);
        Isi = 0.09 * d * f * (Vm - Vsi);

        alpha_d = 0.095*exp(-0.01*(Vm-5))/(1+exp(-0.072*(Vm-5)));
        beta_d = 0.07*exp(-0.017*(Vm+44))/(1+exp(0.05*(Vm+44)));
        alpha_f = 0.012*exp(-0.008*(Vm+28))/(1+exp(0.15*(Vm+28)));
        beta_f = 0.0065*exp(-0.02*(Vm+30))/(1+exp(-0.2*(Vm+30)));

        Ca_conc_dot = -1e-4 * Isi + 0.07 * (1e-4 - Ca_conc);

        % Time-dependent potassium current
        gKmax = 0.282 * sqrt(K_conc / 5.4);
        if (Vm > -100)
            Xi = 2.837 * (exp(0.04*(Vm+77))-1)/((Vm+77)*exp(0.04*(Vm+35)));
        else
            Xi = 1;
        end

        IK = gKmax * X * Xi * (Vm - VK);

        alpha_X = 0.0005 * exp(0.083*(Vm+50))/(1+exp(0.057*(Vm+50)));
        beta_X = 0.0013 * exp(-0.06*(Vm+20))/(1+exp(-0.04*(Vm+20)));

        % Time-independent potassium current
        gKimax = 0.6047 * sqrt(K_conc / 5.4);
        alpha_Ki = 1.02/(1+exp(0.2385*(Vm-VKi-59.215)));
        beta_Ki = (0.49124*exp(0.08032*(Vm-VKi+5.476))+exp(0.06175*(Vm-VKi-594.31)))/(1+exp(-0.5143*(Vm-VKi+4.753)));
        Ki_inf = alpha_Ki / (alpha_Ki + beta_Ki); 
        IKi = gKimax * Ki_inf * (Vm - VKi);

        % Plateu potassium current
        VKp = VKi;
        Kp = 1/(1+exp((7.488-Vm)/5.98));
        IKp = 0.0183 * Kp * (Vm - VKp);
        
        % Background current
        Ib = 0.03921*(Vm+59.87);

        I_tot = Ib + IKp + IKi + IK + Isi + INa;

        Vdot = (1/Cm) * (I_t(t) - I_tot); % V/s
        hdot = alpha_h * (1 - h) - beta_h * h;
        jjdot = alpha_jj * (1 - jj) - beta_jj * jj;
        mdot = alpha_m * (1 - m) - beta_m * m;
        ddot = alpha_d * (1 - d) - beta_d * d;
        fdot = alpha_f * (1 - f) - beta_f * f;
        Xdot = alpha_X * (1 - X) - beta_X * f;
        Kidot = alpha_Ki * (1 - Ki) - beta_Ki * Ki;
        XXdot = [Vdot; Ca_conc_dot; hdot; jjdot; mdot; ddot; fdot; Xdot; Kidot];
%         if (~isreal(XXdot))
%             pause;
%         end
%         pause;

    end

    
end