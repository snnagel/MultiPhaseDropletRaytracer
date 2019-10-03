function [ ] = drawPlane( p, sl )
% Draws a plane p with equation p(2)*x1+p(3)*x2+p(4)*x3= p(1)
%centers around point closest to origin and draws side length scaled by sl

n=p(2:end); %normal


% 2D
if length(p)==3
    p1=p(1)*n+[-n(2), n(1)]*sl/2
    p2=p(1)*n-[-n(2), n(1)]*sl/2
    plot([p1(1), p2(1)], [p1(2), p2(2)], 'b')
end

% 3D
if length(p)==4
x = [1 -1 -1 1]*sl; % Generate data for x vertices
y = [1 1 -1 -1]*sl; % Generate data for y vertices
z = -1/p(4)*(p(2)*x + p(3)*y - p(1)) % Solve for z vertices data
patch('XData',x,'YData',y,'ZData',z)
end




end

