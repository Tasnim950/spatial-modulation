function [omega] = omega36(xr, xi, Nr)
m = length(xr);
s = 0;

for i1 = 1:m
    for i1p = 1:m
        for i2 = 1:m
            for i2p = 1:m
                if (i1 == i1p && i2 == i2p)
                    continue;
                end
                s1 = (xr(i1) - xr(i1p))^2 + (xi(i1) - xi(i1p))^2;
                s2 = (xr(i2) - xr(i2p))^2 + (xi(i2) - xi(i2p))^2;
                
                s = s + (s1 + s2)^(-Nr);
            end
        end
    end
end
omega = s;
end