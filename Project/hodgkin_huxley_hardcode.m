function II = hodgkin_huxley_hardcode(tf, I_t, params)
% Source: these parameters come from https://github.com/swharden/pyHH

%   For this model, the time scale is ms

    Cm = 10 * 1e-3; % nF/mm^2 -> uF/mm^2
    gKmax = 50; % uS/mm^2
    gNamax = 1000; % uS/mm^2 
    gL = 3; % uS/mm^2
    VK = -35; % mV
    VNa = 115; % mV
    VL = 10.6; % mV

%     Cm = 1 * 1e-3; % uF/cm^2
%     gKmax = 35; % mS/cm^2
%     gNamax = 40; % mS/cm^2 
%     gL = 0.3; % mS/cm^2
%     VK = -77; % mV
%     VNa = 55; % mV
%     VL = -65; % mV
    opt = odeset('RelTol', 1e-3, 'MaxStep', 0.05);
    initX = fsolve(@(x) simX(0, x), [0;0.5;0.5;0.5]);
%     initX = [-65; 0; 0; 0];
    initR = getR(initX(1));
    % The initial values of the m, n, and h
    initU = [initR(1)/(initR(1)+initR(2));  ...
                initR(3)/(initR(3)+initR(4));  ...
                initR(5)/(initR(5)+initR(6))];

    [tt, XX] = ode45(@simX, [0, 100], [0; initU], opt);
    figure(1)
    plot(tt, XX(:, 1), "b-")
    title(sprintf("Cm = %.2f", Cm), "Interpreter", "latex")
    figure(2)
    plot(tt, XX(:, 2), "r-")
    hold on;
    plot(tt, XX(:, 3), "g-")
    plot(tt, XX(:, 4), "m-")
    legend("n", "m", "h")
    title(sprintf("Cm = %.2f", Cm), "Interpreter", "latex")
    hold off;

    function Xdot = simX(t, X)
        Vm = X(1); % mV
        U = X(2:4);
        I_K = gKmax * (U(1) .^ 4) * (Vm - VK); % mA/mm^2 = S/mm^2 * dimless * mV
        I_Na = gNamax * (U(2) .^ 3) * U(3) * (Vm - VNa);
        I_L = gL * (Vm - VL);
        I_ion = I_K + I_Na + I_L;
        Vmdot = (I_t(t) * 1e6 - I_ion) / Cm; % nA/mm^2 / F/mm^2 = mV/s
        % mV/s -> mV/ms: / 1000
 
        R = getR(Vm);
        Udot = R(1:2:5) .* (1 - U) - R(2:2:6) .* U;
        % Udot has the unit of /ms

        Xdot = [Vmdot; Udot];
    end

    function R = getR(Vm) 
        % the parameter A is in the unit of /ms
        alpha_n = .01 * ((10-Vm) / (exp((10-Vm)/10)-1));
        beta_n = .125*exp(-Vm/80);
        alpha_m = .1*((25-Vm) / (exp((25-Vm)/10)-1));
        beta_m = 4*exp(-Vm/18);
        alpha_h = .07*exp(-Vm/20);
        beta_h = 1/(exp((30-Vm)/10)+1);
%         Vm = Vm + 1e-8;
%         alpha_n = 0.02*(Vm-25)/(1-exp(-(Vm-25)/9));
%         beta_n = -0.002*(Vm-25)/(1-exp((Vm-25)/9));
%         alpha_m = 0.182*(Vm+35)/(1-exp(-(Vm+35)/9));
%         beta_m = -0.124*(Vm+35)/(1-exp((Vm+35)/9));
%         alpha_h = 0.25*exp(-(Vm+90)/12);
%         beta_h = 0.25*exp((Vm+62)/6)/exp((Vm+90)/12);       
        
        R = [alpha_n; beta_n; alpha_m; beta_m; alpha_h; beta_h]; % 1/ms
    end
    
end