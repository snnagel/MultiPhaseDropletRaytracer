%3D version of Generate Ray paths- doesn't include coverslip refraction

% Make sure DropletTracer folder is in your path. 
addpath('DropletTracer')

%Droplet parameters:
vr=1; %volume ratio
Rd=50; %droplet radius
Ri=-50; %internal radius of curvature
dropLocation=[0,0,0];
dropAngle=0;  % rotation angle (degrees)

% Indexes of refraction: 
nm=1.33; % surrounding Medium
no=1.27;% innerPhase
ni=1.387; % outerPhase

%set number of rays: 
rayDensity=0.05 %rays per unit length

%set ray start points: 
N=round(2*Rd*rayDensity)% total number of rays is N^2
xStart=linspace(-Rd, Rd, N) %full width of droplet,
yStart=linspace(-Rd, Rd, N) %full width of droplet,
numRays=length(xStart);
zStart=(Rd+5); %5 units above droplet. 

%set ray start directions: 
rDir=[0, 0, -1];

%Ray amplitude is between 0 and 1, through away rays that have less than:
minimumAmplitude=0.1;

%plane to propagate rays to after interaction with the droplet
%if rays don't hit this plane they are thrown away:
%equation of plane is p(1)=p(2)x+p(3)y. 
screen=[ -300,0,0, 1]; % Plane that will be the screen.   y=-300
screenSize=100
%you can choose to not use a screen: 
useScreen=true; 
propagationDistance=100; %instead propagate all rays after they leave the drop by this amount

n=1;

% From contact angle and volume ratio, determine radius of curvature: 
    
    % Make the drop object:
    theDrop=drop(dropLocation, Rd, Ri, no, ni, vr, dropAngle);

    % An Array to store the Rays in: 
    Rays=[]; %This is slow, but for this application its ok.

    % Send Rays through scene
    for ii=1:N
        for jj=1:N
        %create ray object: 
        if xStart(ii)^2+yStart(jj)^2<Rd^2 %only those hitting drop
        aRay=ray([xStart(ii), yStart(jj), zStart], rDir, nm);
        
        %this is the line that actually does the raytracing calculation:
        aRayTraced=traceRay(aRay, theDrop, nm, minimumAmplitude);
        
        %propagate to screen and store
        for t=aRayTraced
            if useScreen
            [loc, norm]=rayPlaneIntersection(t, screen);
            if (abs(loc(1))< screenSize)
                t=propagate(t, loc);
                Rays=[Rays, t];
            end
            else
                t=propagate(t, propagationDistance);
                Rays=[Rays, t];
            end
            
       end
        end
        end
    end
    
    %% plot the output: 
    figure
    % Draw the rays.  (Modify to draw only selected rays.)
    for Ray=Rays
        drawRay(Ray);
        hold on
    end
       
    drawDrop(theDrop) 
    %if you wish to modify colors of the drop, modify drawDrop() in Drop.m  
    
    axis image
    xlabel('x')
    ylabel('y')
    zlabel('z')
    %axis off