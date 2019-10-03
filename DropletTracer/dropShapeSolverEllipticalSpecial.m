function y=dropShapeSolverEllipticalSpecial(x, Rd, a, b, vr)
d=x(1);
L=x(2);
Ha=Rd-L;
Hb=d-Rd+L;
Hc=Rd+L;
Vh=Ha^3/3*(3*Rd-Ha)+b^2*(Hb^3/(3*a^2)-Hb+2/3*a);
Vf=Hc^3/3*(3*Rd-Hc)-b^2*(Hb^3/(3*a^2)-Hb+2/3*a);

y(3)=Vh-(vr/(vr+1))^(1/3);
y(2)=Vf-(1/(vr+1))^(1/3);
y(1)=1-(L-d)^2/a^2-(Rd^2-L^2)/b^2;
%y(1)=(L*(b^2 - a^2)-b^2*d)^2-((-Rd^2 + d^2)*a^2*b^2 + a^2*b^4 + Rd^2*a^4 - a^4*b^2);
%y(2)=2*a^3*b^2*(1 + vr) + b^2*(d - L)^3*(1 + vr) + a^2*(-2*Rd^3*(-1 + vr) - 3*b^2*(d - L)*(1 + vr) + L^3*(1 + vr) - 3*L*Rd^2*(1 + vr));
end