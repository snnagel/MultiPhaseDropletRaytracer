%% Droplet parameters:
vr=1; %volume ratio.  V_internal/v_external
Rd=50; %droplet radius
dropLocation=[100,50]%, 0];
dropAngle=0;%45;  % rotation angle (degrees)
Ri1=50; %Internal radius of curvature. use positve numbers
Ri2=50

% Indexes of refraction: 
nm=1.33; % Surrounding Medium
n1=1.4
n2=1.27; 
n3=1.39;
v1=0.3
v2=0.3

theDrop=TriPhaseDrop(dropLocation, Rd, Ri1, Ri2, n1, n2, n3, v1, v2, dropAngle);
figure
drawDrop(theDrop) 
axis image
set(gcf, 'color', 'white')
axis off
