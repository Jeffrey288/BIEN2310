function [] = stirredtank_Jeffrey_HartantoKwee(Cin, tau, tf, h)
    function val = stirredtank_func(~, C)
        val = 1/tau * (Cin - C);
    end
    [ftt, fyy] = feuler(@stirredtank_func, 0, tf, h);
    [btt, byy] = beuler_Jeffrey_HartantoKwee(@stirredtank_func, 0, tf, h);
    figure;
    plot(ftt, fyy, "b-", btt, byy, "r-");
    title('HARTANTO KWEE, Jeffrey. HW2 Q3(d)')
end


%{
For part (d),
The backward Euler seems to be more acccurate.
With the following parameters, it is observed that the foward Euler
solution behavies erratically:
    stirredtank_Jeffrey_HartantoKwee(1, 0.5, 5, 0.8)
The backward Euler behaves well and stable because solving the implicit equation exactly using fzero() helps us find an exactly value such that dy/dt matches the prediction by y_{i+1}=y_i+\Delta t*f(t,y). Solving this equation exactly enables us to use a larger step size without worrying about stability. (The reason we use a small step size is because the linear approximation is based on the assumption that \Delta t is small. The found  value every step is very close to the real value, so the error accumulated at each step is small and prevents the apporximator from fluctuating widly.
%}