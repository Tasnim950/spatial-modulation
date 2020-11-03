function [omega] = omega21(xr, xi, Nr)
m = length(xr);
s = 0;

for i1 = 1:m
    for i1p = 1:m
        for i2 = 1:m
            for i2p = 1:m
                s = s + (xr(i1p)^2 + xr(i2)^2 + (xr(i2p) - xr(i1))^2 ...
                    + xi(i1)^2 + xi(i1p)^2 + (xi(i2) - xi(i2p))^2)^(-Nr);
            end
        end
    end
end
omega = s; 
end