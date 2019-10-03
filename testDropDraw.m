%% Droplet parameters:
vr=1; %volume ratio.  V_internal/v_external
Rd=50; %droplet radius
dropLocation=[100,50]%, 0];
dropAngle=0;45;  % rotation angle (degrees)
Ri=50; %Internal radius of curvature. use positve numbers

% Indexes of refraction: 
nm=1.33; % Surrounding Medium
no=1.27; % Outer Phase
ni=1.39;% Inner Phase

theDrop=drop(dropLocation, Rd, Ri, no, ni, vr, dropAngle);
figure
drawDrop(theDrop) 
axis image
set(gcf, 'color', 'white')
axis off
