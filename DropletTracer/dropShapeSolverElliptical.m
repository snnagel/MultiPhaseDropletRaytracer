function y=dropShapeSolverElliptical(x, Rd, a, b, vr)
d=x(1);
L=x(2);
y(1)=(L*(b^2 - a^2)-b^2*d)^2-((-Rd^2 + d^2)*a^2*b^2 + a^2*b^4 + Rd^2*a^4 - a^4*b^2);
y(2)=2*a^3*b^2*(1 + vr) + b^2*(d - L)^3*(1 + vr) + a^2*(-2*Rd^3*(-1 + vr) - 3*b^2*(d - L)*(1 + vr) + L^3*(1 + vr) - 3*L*Rd^2*(1 + vr));
end