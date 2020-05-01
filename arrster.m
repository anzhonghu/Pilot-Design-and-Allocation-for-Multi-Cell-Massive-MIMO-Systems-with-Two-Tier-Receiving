function [a] = arrster(x, y, M, Nh, Nv, BS_height)
a = zeros(M, 1);
phi = atan(y/x);
if phi<=0
    phi = pi/3 + phi;
else
    phi = -pi/6-(pi/2-phi);
end
theta = pi/2 - atan(sqrt(x^2+y^2)/BS_height);
for m = 1 : Nh
    for n = 1 : Nv
        a((m-1)*Nv+n, 1) = 1 / sqrt(M) * exp(1i*pi*(-(n-1)*sin(theta)+sin(phi)*cos(theta)*(m-1)));
    end
end

end

