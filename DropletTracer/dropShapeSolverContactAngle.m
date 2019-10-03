function y=dropShapeSolverContactAngle(x, Rd, vr, contactAngle)
    %contact angle is the COSINE OF THE CONTACT ANGLE

    x1=x(1);
    Ri=x(2);
    y1 = (x1+Rd-Ri).^2.*(x1.^2+2.*x1.*(Ri-Rd)-3.*(Rd+Ri).^2)+16*x1.*Rd.^3./(1+vr);
    y2= contactAngle+(x1^2-Rd^2-Ri^2)/(2*Ri*Rd);
    y=[y1, y2];
end