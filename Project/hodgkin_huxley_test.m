% [Cm, gKmax, gNamax, gL, VK, VNa, VL] = params(1:7);
% Uparam = params(8:(8+24-1));

% hodgkin_huxley(30, @I_none, ...
% [...
%     [10; 360; 1200; 3; -77; 50; -54.402]; ...
%     [...
%         0.1 * 40; 0.1; 1; 40; -1/0.1; -1;...
%         4; 0; 0; 65; 1/0.0556; 1;...
%         0.01 * 55; 0.01; 1; 55; -1/0.1; -1;...
%         0.125; 0; 0; 65; 1/0.0125; 1;...
%         0.07; 0; 0; 65; 1/0.05; 1;...
%         1; 0; 1; 35; -1/0.1; 1;...
%     ]...
% ])

% hodgkin_huxley_hardcode_voltage_clamp(100, @(x) -65 + (x > 20) * (65 - 45), NaN)

% Parameters adapted from
% https://neuronaldynamics.epfl.ch/online/Ch2.S2.html
% These parameters are for "pyramidal neurons of the cortex"
% mS/cm^2 = 10 uS/mm^2
% hodgkin_huxley(100, @I_step, ...
% [...
%     [10; 350; 400; 30; -77; 55; -65]; ...
%     [...
%         0.02 * -25; 0.02; 1; -25; -9; -1;...
%         -0.002 * -25; -0.002; 1; -25; 9; -1;...
%         0.182 * 35; 0.182; 1; 35; -9; -1;...
%         -0.124 * 35; -0.124; 1; 35; 9; -1;...
%         0.25; 0; 0; 90; 12; 1;...
%         0.25 * exp(17/6); 0; 0; 0; -12; 1;...
%     ]...
% ], "Vm0", -66);
% , "Vm0", -61

% hodgkin_huxley_hardcode_normal(300, @I_step, [0, 1, 2], "Vm0", -65);


% function I = I_const(t)
%     I = 0; % 2.7 uA/cm^2 -> 2.7e-5 mA/cm^2
% end
% 
% function I = I_step(t)
%     if (t >= 20)
%         I = 0.31e-6; % aspecify in terms of mA/mm^2
%     else 
%         I = 0;
%     end
% end

I_i = 7.47; % in uA/cm^2
accuracy = 0; % try I_i = 0.37 with Fast and you'll know why I added this
tf = 600; % in ms
I_step = @(t) ((t > 20) * (t < 180)) * 15; % the injection current has the unit uA/cm^2
% hodgkin_huxley_hardcode_original(tf, @(t) (t > 20) * 7.45, NaN, "Vm0", 0, 'accuracy', 0);
% hodgkin_huxley_hardcode_original(tf, I_step, NaN, "Vm0", 0, 'accuracy', accuracy);
hodgkin_huxley_hardcode_LRd(tf, I_step, NaN, 'accuracy', accuracy);

% N = 50; % number of points to calculate
% Imax = 20; % range of amplitude to try
% II = linspace(0, Imax, N);
% freq = zeros(N, 1);
% for i = 1:N
%     I_step = @(t) (t > 20) * II(i);
%     freq(i) = hodgkin_huxley_hardcode_original(300, I_step, NaN, "Vm0", 0, 'accuracy', 1, 'getFreq', true);
% end
% figure(1);
% plot(II, freq, "b-");
% title("Spike Frequency vs Injected Current")
% xlabel("injected current (mA/mm^2)")
% ylabel("frequency (Hz)")
% ylim([0, 100])

% https://neuronaldynamics.epfl.ch/online/Ch2.S2.html
% https://bernstein-network.de/wp-content/uploads/2021/02/04_Lecture-04-Hodgkin-Huxley-model.pdf