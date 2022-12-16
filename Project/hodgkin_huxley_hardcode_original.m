function res = hodgkin_huxley_hardcode_original(tf, I_t, params, options)

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

    Cm = 1; % uF/cm^2
    gKmax = 36; % mS/cm^2
    gNamax = 120; % mS/cm^2 
    gL = 0.3; % mS/cm^2
    VK = -12; % mV
    VNa = 115; % mV
    VL = 10.613; % mV

    if isnan(options.Vm0)
        initX = fsolve(@(x) simX(0, x), [0; 0.5; 0.5; 0.5]);
    else
%         % The initial values of the m, n, and h
        initR = getR(options.Vm0);
        initU = [ ...
                    initR(1)/(initR(1)+initR(2));  ...
                    initR(3)/(initR(3)+initR(4));  ...
                    initR(5)/(initR(5)+initR(6))
                ];
%         initU(3) = 1/(1+exp(options.Vm0-(-65)*6.2));
        initX = [options.Vm0; initU];
    end

    if (options.accuracy == 0)
        opt = odeset('RelTol', 1e-3, 'MaxStep', 1.0);
    elseif (options.accuracy == 1)
        opt = odeset('RelTol', 1e-5, 'MaxStep', 0.2);
    else
        opt = odeset('RelTol', 1e-7, 'MaxStep', 0.1);
    end
    [tt, XX] = ode45(@simX, [0, tf], initX, opt);

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
%     return;
%     figure(3)
%     plot(tt, XX(:, 2), "r-")
%     hold on;
%     plot(tt, XX(:, 3), "g-")
%     plot(tt, XX(:, 4), "m-")
%     legend("$n$", "$m$", "$h$", "Interpreter", "latex")
%     title("Subunit Activations", "Interpreter", "latex")
%     xlabel("time (ms)")
%     ylabel("subunit activation")
%     hold off;
%     figure(4)
%     I_K = gKmax * (XX(:, 2) .^ 4) .* (XX(:, 1) - VK);
%     I_Na = gNamax * (XX(:, 3) .^ 3) .* XX(:, 4) .* (XX(:, 1) - VNa);
%     subplot(2,1,1);
%     plot(tt, I_K, "g-")
%     subplot(2,1,2);
%     title("Potassium Ion Current", "Interpreter", "latex")
%     xlabel("time (ms)")
%     ylabel("Current (A/cm^2)")
%     legend("$I_K$", "Interpreter", "latex")
%     plot(tt, I_Na, "m-")
%     legend("$I_{Na}$", "Interpreter", "latex")
%     title("Sodium Ion Current", "Interpreter", "latex")
%     xlabel("time (ms)")
%     ylabel("Current (A/cm^2)")
%     hold off;
%     
    function Xdot = simX(t, X)

        Vm = X(1); % mV
        U = X(2:4);
        I_K = gKmax * (U(1) .^ 4) * (Vm - VK); % uA/cm^2 = mS/cm^2 * dimless * mV
        I_Na = gNamax * (U(2) .^ 3) * U(3) * (Vm - VNa);
        I_L = gL * (Vm - VL);
        I_ion = I_K + I_Na + I_L;
        Vmdot = (I_t(t) - I_ion) / Cm; % uA/cm^2 / uF/cm^2 = V/s -> mV/ms
 
        R = getR(Vm);
        Udot = R(1:2:5) .* (1 - U) - R(2:2:6) .* U;
        % Udot has the unit of /ms

        Xdot = [Vmdot; Udot];
    end

    function R = getR(Vm) 
        % the parameter A is in the unit of /ms
        Vm = Vm + 1e-8; % prevent division by 0
        alpha_n = 0.01*(-Vm+10)/(exp((-Vm+10)/10)-1);
        beta_n = 0.125*exp(-Vm/80);
        alpha_m = 0.1*(-Vm+25)/(exp((-Vm+25)/10)-1);
        beta_m = 4*exp(-Vm/18);
        alpha_h = 0.07*exp(-Vm/20);
        beta_h = 1/(exp((-Vm+30)/10)+1);
        R = [alpha_n; beta_n; alpha_m; beta_m; alpha_h; beta_h]; % 1/ms
    end
    
end