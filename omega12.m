function [omega] = omega12(xr, xi, Nr)
m = length(xr);
s = 0; 

for i1 = 1:m
    for i1p = 1:m
        s1 = xr(i1)^2 + xr(i1p)^2 + (xi(i1) - xi(i1p))^2;
        for i2 = 1:m
            for i2p = 1:m
                s2 = (xr(i2) - xr(i2p))^2 + xi(i2)^2 + xi(i2p)^2;
                
                s = s + (s1 + s2)^(-Nr);
            end
        end
    end
end
omega = s;
end