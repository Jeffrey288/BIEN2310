function II = hodgkin_huxley_hardcode_voltage_clamp(tf, V_t, params, options)

    arguments
        tf
        V_t
        params
        options.step (1,1) int8 = true % if the function is a step function, set this as yes
        options.accuracy (1,1) int8 = 0
    end

%   For this model, the time scale is ms

    Cm = 1 * 1e-3; % uF/cm^2
    gKmax = 35; % mS/cm^2
    gNamax = 40; % mS/cm^2 
    gL = 0.3; % mS/cm^2
    VK = -77; % mV
    VNa = 55; % mV
    VL = -65; % mV

    % The initial values of the m, n, and h
    initR = getR(V_t(0));
    initU = [ ...
                initR(1)/(initR(1)+initR(2));  ...
                initR(3)/(initR(3)+initR(4));  ...
                initR(5)/(initR(5)+initR(6))
            ];


    if (options.accuracy == 0)
        opt = odeset('RelTol', 1e-3, 'MaxStep', 0.5);
    elseif (options.accuracy == 1)
        opt = odeset('RelTol', 1e-5, 'MaxStep', 0.2);
    else
        opt = odeset('RelTol', 1e-7, 'MaxStep', 0.1);
    end

    [tt, UU] = ode45(@simX, [0, tf], initU, opt);

    I_K = gKmax * (UU(:, 1) .^ 4) .* (V_t(tt) - VK); % nA/mm^2 = uS/mm^2 * dimless * mV
    I_Na = gNamax * (UU(:, 2) .^ 3) .* UU(:, 3) .* (V_t(tt) - VNa);
    I_L = gL * (V_t(tt) - VL);
    I_ion = I_K + I_Na + I_L;
    I_t = (dVdt(tt) * Cm * 1000 + I_ion); % mA/mm^2 / F/mm^2 = V/s

    hold off;
    figure(1)
    plot(tt, arrayfun(V_t, tt), "c-")
    title("Applied Voltage", "Interpreter", "latex")
    xlabel("time (ms)")
    ylabel("voltage (mV)")

    figure(2)
    plot(tt, I_t, "b-");
    title("Resultant Current", "Interpreter", "latex")
    xlabel("time (ms)")
    ylabel("current (mA/mm^2)")

    figure(3)
    plot(tt, UU(:, 1), "r-")
    hold on;
    plot(tt, UU(:, 2), "g-")
    plot(tt, UU(:, 3), "m-")
    legend("$n$", "$m$", "$h$", "Interpreter", "latex")
    title("Subunit Activations", "Interpreter", "latex")
    xlabel("time (ms)")
    ylabel("subunit activation")
    hold off;

    figure(4)
    I_K = gKmax * (UU(:, 1) .^ 4) .* (V_t(tt) - VK);
    I_Na = gNamax * (UU(:, 2) .^ 3) .* UU(:, 3) .* (V_t(tt) - VNa);
    subplot(2,1,1);
    plot(tt, I_K, "g-")
    subplot(2,1,2);
    title("Potassium Ion Current", "Interpreter", "latex")
    xlabel("time (ms)")
    ylabel("Current (nA/mm^2)")
    legend("$I_K$", "Interpreter", "latex")
    plot(tt, I_Na, "m-")
    legend("$I_{Na}$", "Interpreter", "latex")
    title("Sodium Ion Current", "Interpreter", "latex")
    xlabel("time (ms)")
    ylabel("Current (nA/mm^2)")
    hold off;
    
    function res = dVdt(t)
        if (options.step)
            res = 0;
        else
            % https://en.wikipedia.org/wiki/Finite_difference_coefficient
            res = 1/12 * V_t(t + -2e-3) - 2/3 * V_t(t + -1e-3) + 2/3 * V_t(t + 1e-3) - 1/12 * V_t(t + 2e-3);
            res = res / 1e-3;
        end
    end

    function Udot = simX(t, U)
        
        Vm = V_t(t);
        R = getR(Vm);
        Udot = R(1:2:5) .* (1 - U) - R(2:2:6) .* U;
        % Udot has the unit of /ms

    end

    function R = getR(Vm) 
        % the parameter A is in the unit of /ms
        Vm = Vm + 1e-8; % prevent division by 0
        alpha_n = 0.02*(Vm-20)/(1-exp(-(Vm-20)/9));
        beta_n = -0.002*(Vm-20)/(1-exp((Vm-20)/9));
        alpha_m = 0.182*(Vm+35)/(1-exp(-(Vm+35)/9));
        beta_m = -0.124*(Vm+35)/(1-exp((Vm+35)/9)); 
        alpha_h = 0.25*exp(-(Vm+90)/12);
        beta_h = 0.25*exp((Vm+62)/6)/exp((Vm+90)/12);       
        R = [alpha_n; beta_n; alpha_m; beta_m; alpha_h; beta_h]; % 1/ms
    end
    
end