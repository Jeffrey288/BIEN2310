function [] = fiboplot_Jeffrey_HartantoKwee()
    F = fiboseries(20);
    ratios = F(2:20, 1) ./ F(1:19, 1);
    args = 1:1:19
    figure;
    plot(args, ratios, 'bx')
    title('HARTANTO KWEE, Jeffrey. HW2 Q1(b)')
end