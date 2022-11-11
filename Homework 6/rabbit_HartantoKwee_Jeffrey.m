function N = rabbit_HartantoKwee_Jeffrey(m, p, d, showplot)

% rabbit_HartantoKwee_Jeffrey(60, 1/30, 800, true)

new_rabbit_pairs = zeros(d, 1);
new_rabbit_pairs(1) = 1;

for i=2:d
    for j=1:(i-m)
        for k=1:new_rabbit_pairs(j)
            if rand(1,1) < p
                new_rabbit_pairs(i) = new_rabbit_pairs(i) + 1;
            end
        end
    end
end

total_rabbit_pairs = new_rabbit_pairs;
for i=2:d
    total_rabbit_pairs(i) = total_rabbit_pairs(i) + total_rabbit_pairs(i-1);
end

N = total_rabbit_pairs(d);

if (showplot)
    plot(1:d, total_rabbit_pairs, 'b-');
end

end