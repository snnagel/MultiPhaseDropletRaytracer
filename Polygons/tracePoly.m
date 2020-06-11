function r = tracePoly( r, N, D, n2 )
%trace ray r through 2D polygon with N sides with diameter D (note half polygons so N=3 is
%a hexagon. 
%this is a quickly done function
    %assume r.y>0 means ray is above polygon
    %r.y< 0 in polygon
    
tol=1e-14 %consider numbers this close to be equal

%angle between planes: 
p_angle=pi/N;
ROT=[cos(p_angle/2), -sin(p_angle/2); sin(p_angle/2), cos(p_angle/2)]
%distance of each side from center: 
c=D/2*cos(p_angle/2);
%half of each side:
b=D/2*sin(p_angle/2);

%is ray in polygon?  if not will it hit polygon?p(2)*x1+p(3)*x2+p(4)*x3= p(1)
   topPlane=[0, 0, 1];
if (r.location(2)>0)
   [loc, norm]=rayPlaneIntersection(r, topPlane);
   if abs(loc(1))< D/2
       r=propagate(r, loc)
   else
       'Ray Missed'
       return
   end   
end

%for each side does the ray hit the side? if so propagate to side and
%reflect
intFound=false;
for ii= 1:N
%  ii=ii
%  rotMat=ROT^(2*(ii-1)+1)
    inorm=ROT^(2*(ii-1)+1)*[1;0]
    iplane=[-c, inorm']
    iAngle=p_angle/2*(2*(ii-1)+1)+pi;
    icenter=[c*cos(iAngle), c*sin(iAngle)];
 %   plot(icenter(1), icenter(2), 'or')
    [loc, norm]=rayPlaneIntersection(r, iplane)
% 'Which Not Workin'
% 
%     not(isnan(loc(1)))
%     sum((loc-icenter).^2)<b^2
% abs(sum(loc-r.location))>tol
% r.location
% abs(sum(loc-r.location))
% tol

    if( not(isnan(loc(1))) && sum((loc-icenter).^2)<b^2 && abs(sum((loc-r.location).^2))>tol)
        r=propagate(r, loc)
        [~, r]=refract(r, norm, n2);
        intFound=true;
        break
    end
end
if (not(intFound))
    message='DIDN"T FIND INTERSECTION'
    return
end

%Does it exit?  If not recursive. 
[loc, norm]=rayPlaneIntersection(r, topPlane);
if abs(loc(1))< D/2
     r=propagate(r, loc);
 else
     r=tracePoly( r, N, D, n2 );
end
end