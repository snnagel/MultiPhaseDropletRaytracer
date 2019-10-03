function [ ] = DrawPoly( N, D )
%Draw bottom half of polygon with 2N sides (draw N sides)
p_angle=pi/N;
theta=[pi:p_angle:2*pi]

hold on
plot(D/2*cos(theta), D/2*sin(theta), 'k', 'linewidth', 2)
axis equal

end

