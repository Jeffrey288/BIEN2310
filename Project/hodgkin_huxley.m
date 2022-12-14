function [] = hodgkin_huxley(tf, I_t, params, options)
    
    arguments
        tf
        I_t
        params
        options.Vm0 (1,1) double = NaN
    end

    % For this model, the time scale is ms

    % Parameters:
    % Cm (nF/mm^2): membrane capacitance
    % gKmax (uS/mm^2): max potassium ion conductance
    % gNamax (uS/mm^2): max sodium ion conductance
    % gL (uS/mm^2): max leakage conductance
    % VK (mV): potassium reversal voltage
    % VNa (mV): sodium reversal voltage
    % VL (mV): leakage potential
    % Uparam: [A (/ms), B, C, D, F, H] for each of R's elements

    % t is in ms
    % Vm is in mV
    % I_t is given in mA/mm^2

%   [Cm, gKmax, gNamax, gL, VK, VNa, VL] = params(1:7)
    Cm = params(1) * 1e-3; % nF/mm^2 -> uF/mm^2
    gKmax = params(2); % uS/mm^2
    gNamax = params(3); % uS/mm^2 
    gL = params(4); % uS/mm^2
    VK = params(5); % mV
    VNa = params(6); % mV
    VL = params(7); % mV
    Uparam = params(8:(8+36-1));
    Uparam = reshape(Uparam, [6, 6]).';
    
    % X stores [Vm; U]

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
        initX = [options.Vm0; initU];
    end
    
    opt = odeset('RelTol', 1e-7, 'MaxStep', 5e-3)
    [tt, XX] = ode45(@simX, [0, tf], initX, opt);

    figure(1)
    plot(tt, XX(:, 1), "b-")
    title("Hodgkin Huxley Model: Membrane Voltage", "Interpreter", "latex")
    xlabel("time (ms)")
    ylabel("membrane potential difference (mV)")
    figure(2)
    plot(tt, XX(:, 2), "r-")
    hold on;
    plot(tt, XX(:, 3), "g-")
    plot(tt, XX(:, 4), "m-")
    legend("$n$", "$m$", "$h$", "Interpreter", "latex")
    title("Subunit Activations", "Interpreter", "latex")
    xlabel("time (ms)")
    ylabel("subunit activation")
    hold off;

    function Xdot = simX(t, X)
        Vm = X(1); % mV
        U = X(2:4);
        I_K = gKmax * (U(1) .^ 4) * (Vm - VK); % nA/mm^2 = uS/mm^2 * dimless * mV
        I_Na = gNamax * (U(2) .^ 3) * U(3) * (Vm - VNa);
        I_L = gL * (Vm - VL);
        I_ion = I_K + I_Na + I_L;
        Vmdot = (I_t(t) * 1e6 - I_ion) / Cm / 1000; % nA/mm^2 / uF/mm^2 = mV/s
        % note that I_t(t) is given as mA/mm^2
        % V/s -> V/ms: / 1000
 
        R = getR(Vm);
        Udot = R(1:2:5) .* (1 - U) - R(2:2:6) .* U;
        % Udot has the unit of /ms

        Xdot = [Vmdot; Udot];
    end

    function R = getR(Vm)         
%         alpha_n = 0.02*(Vm-25)/(1-exp(-(Vm-25)/9));
%         beta_n = -0.002*(Vm-25)/(1-exp((Vm-25)/9));
%         alpha_m = 0.182*(Vm+35)/(1-exp(-(Vm+35)/9));
%         beta_m = -0.124*(Vm+35)/(1-exp((Vm+35)/9));
%         alpha_h = 0.25*exp(-(Vm+90)/12);
%         beta_h = 0.25*exp((Vm+62)/6)/exp((Vm+90)/12);       
%         R = [alpha_n; beta_n; alpha_m; beta_m; alpha_h; beta_h] % 1/ms
        Vm
        R = (Uparam(:, 1) + Uparam(:, 2) .* (Vm + 1e-8)) ./ ...
            (Uparam(:, 3) + Uparam(:, 6) .* exp(((Vm + 1e-8) + Uparam(:, 4)) ./ Uparam(:, 5)))
    end


end