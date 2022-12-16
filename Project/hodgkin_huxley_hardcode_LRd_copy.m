function res = hodgkin_huxley_hardcode_LRd_copy(tf, I_t, params, options)
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
%     initX = zeros(9,1);
%     initX(2) = 4e-4;
%     initX = fsolve(@(x) simX(0, x), initX);
%     initX(2) = 4e-4;
    initX = [0.982660523699656, 0.989108212766685, 0.00171338077730188, 0.00017948816388306, -84.3801107371, 0.00302126301779861, 0.999967936476325, 0.0417603108167287];

    if (options.accuracy == 0)
        opt = odeset('MaxStep', 0.1);
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
    plot(tt, XX(:, 5), "b-")
    title("Original HH: Membrane Voltage", "Interpreter", "latex")
    subtitle(sprintf("spike period = %.1f ms", period), "Interpreter", "latex")
    xlabel("time (ms)")
    ylabel("membrane potential difference (mV)")
    ylim([-30, 130])

    function dY = simX(t, Y)

        E_b = -59.87;   % millivolt (in background_current)
        g_b = 0.03921;   % milliS_per_cm2 (in background_current)
        g_Na = 23.0;   % milliS_per_cm2 (in fast_sodium_current)
        Ki = 145.0;   % millimolar (in ionic_concentrations)
        Ko = 5.4;   % millimolar (in ionic_concentrations)
        Nai = 18.0;   % millimolar (in ionic_concentrations)
        Nao = 140.0;   % millimolar (in ionic_concentrations)
        C = 1.0;   % microF_per_cm2 (in membrane)
        F = 96484.6;   % coulomb_per_mole (in membrane)
        R = 8314.0;   % joule_per_kilomole_kelvin (in membrane)
        T = 310.0;   % kelvin (in membrane)
        stim_amplitude = -25.5;   % microA_per_cm2 (in membrane)
        stim_duration = 2.0;   % millisecond (in membrane)
        stim_end = 9000000.0;   % millisecond (in membrane)
        stim_period = 1000.0;   % millisecond (in membrane)
        stim_start = 100.0;   % millisecond (in membrane)
        g_Kp = 0.0183;   % milliS_per_cm2 (in plateau_potassium_current)
        PR_NaK = 0.01833;   % dimensionless (in time_dependent_potassium_current)
                
        % time (millisecond)
        
        i_b = g_b*(Y(5)-E_b);
        E_Na = R*T/F*log(Nao/Nai);
        i_Na = g_Na*Y(3)^3.0*Y(1)*Y(2)*(Y(5)-E_Na);
        
        if (Y(5) < -40.0)
           alpha_h = 0.135*exp((80.0+Y(5))/-6.8);
        else
           alpha_h = 0.0;
        end;
        
        if (Y(5) < -40.0)
           beta_h = 3.56*exp(0.079*Y(5))+310000.0*exp(0.35*Y(5));
        else
           beta_h = 1.0/(0.13*(1.0+exp((Y(5)+10.66)/-11.1)));
        end;
        
        dY(1, 1) = alpha_h*(1.0-Y(1))-beta_h*Y(1);
        
        if (Y(5) < -40.0)
           alpha_j = (-127140.0*exp(0.2444*Y(5))-0.00003474*exp(-0.04391*Y(5)))*(Y(5)+37.78)/(1.0+exp(0.311*(Y(5)+79.23)));
        else
           alpha_j = 0.0;
        end;
        
        if (Y(5) < -40.0)
           beta_j = 0.1212*exp(-0.01052*Y(5))/(1.0+exp(-0.1378*(Y(5)+40.14)));
        else
           beta_j = 0.3*exp(-0.0000002535*Y(5))/(1.0+exp(-0.1*(Y(5)+32.0)));
        end;
        
        dY(2, 1) = alpha_j*(1.0-Y(2))-beta_j*Y(2);
        alpha_m = 0.32*(Y(5)+47.13)/(1.0-exp(-0.1*(Y(5)+47.13)));
        beta_m = 0.08*exp(-Y(5)/11.0);
        dY(3, 1) = alpha_m*(1.0-Y(3))-beta_m*Y(3);
        E_si = 7.7-13.0287*log(Y(4)/1.0);
        i_si = 0.09*Y(6)*Y(7)*(Y(5)-E_si);
        dY(4, 1) = -0.0001/1.0*i_si+0.07*(0.0001-Y(4));
        
        g_K = 0.282*sqrt(Ko/5.4);
        
        if (Y(5) > -100.0)
           Xi = 2.837*(exp(0.04*(Y(5)+77.0))-1.0)/((Y(5)+77.0)*exp(0.04*(Y(5)+35.0)));
        else
           Xi = 1.0;
        end;
        
        E_K = R*T/F*log((Ko+PR_NaK*Nao)/(Ki+PR_NaK*Nai));
        i_K = g_K*Y(8)*Xi*(Y(5)-E_K);
        g_K1 = 0.6047*sqrt(Ko/5.4);
        E_K1 = R*T/F*log(Ko/Ki);
        alpha_K1 = 1.02/(1.0+exp(0.2385*(Y(5)-E_K1-59.215)));
        beta_K1 = (0.49124*exp(0.08032*(Y(5)+5.476-E_K1))+1.0*exp(0.06175*(Y(5)-(E_K1+594.31))))/(1.0+exp(-0.5143*(Y(5)-E_K1+4.753)));
        K1_infinity = alpha_K1/(alpha_K1+beta_K1);
        i_K1 = g_K1*K1_infinity*(Y(5)-E_K1);
        Kp = 1.0/(1.0+exp((7.488-Y(5))/5.98));
        E_Kp = E_K1;
        i_Kp = g_Kp*Kp*(Y(5)-E_Kp);
        dY(5, 1) = -1.0/C*(-I_t(t)+i_Na+i_si+i_K+i_K1+i_Kp+i_b);
        alpha_d = 0.095*exp(-0.01*(Y(5)-5.0))/(1.0+exp(-0.072*(Y(5)-5.0)));
        beta_d = 0.07*exp(-0.017*(Y(5)+44.0))/(1.0+exp(0.05*(Y(5)+44.0)));
        dY(6, 1) = alpha_d*(1.0-Y(6))-beta_d*Y(6);
        alpha_f = 0.012*exp(-0.008*(Y(5)+28.0))/(1.0+exp(0.15*(Y(5)+28.0)));
        beta_f = 0.0065*exp(-0.02*(Y(5)+30.0))/(1.0+exp(-0.2*(Y(5)+30.0)));
        dY(7, 1) = alpha_f*(1.0-Y(7))-beta_f*Y(7);
        alpha_X = 0.0005*exp(0.083*(Y(5)+50.0))/(1.0+exp(0.057*(Y(5)+50.0)));
        beta_X = 0.0013*exp(-0.06*(Y(5)+20.0))/(1.0+exp(-0.04*(Y(5)+20.0)));
        dY(8, 1) = alpha_X*(1.0-Y(8))-beta_X*Y(8);

    end

    
end