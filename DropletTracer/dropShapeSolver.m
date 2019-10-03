function y=dropShapeSolver(x, Rd, Ri, vr)
    y = (x+Rd-Ri).^2.*(x.^2+2.*x.*(Ri-Rd)-3.*(Rd+Ri).^2)+16*x.*Rd.^3./(1+vr);
end