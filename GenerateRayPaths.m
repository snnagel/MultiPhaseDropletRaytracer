%There are three primary components of the raytracer: 
%rays, drops, and the traceRay function. 

%The raytracer is untiless, but easy enough to assume all units are in
%micrometers.

% Make sure DropletTracer folder is in your path. 
addpath('DropletTracer')

%% Droplet parameters:
vr=1; %volume ratio.  V_internal/v_external
Rd=50; %droplet radius
dropLocation=[0,0];
dropAngle=0;  % rotation angle (degrees)
Ri=-60; %Internal radius of curvature. use positve numbers

% Indexes of refraction: 
nm=1.33; % Surrounding Medium
no=1.27; % Outer Phase
ni=1.39;% Inner Phase



%% Raytracing Parameters:
%set number of rays: 
rayDensity=0.2 %rays per unit length

%include coverslip (or any glass substrate)
useCoverslip=true;
if useCoverslip
    cThickness=105%thickness of coverslip
    nc=1.52 %refractive index of coverslip
    na=1% refractive index outside of coverslip (air probably)
end
%if you use a thicker substrate, check the location of the screen so that it is after the bottom edge of the coverslip
screenDistance=Rd+1000;
%plane to propagate rays to after interaction with the droplet
%if rays don't hit this plane they are thrown away:
%equation of plane is p(1)=p(2)x+p(3)y. 
screen=[ -screenDistance,0, 1]; % Plane that will be the screen.  
screenSize=200
%you can choose to not use a screen: 
useScreen=true; 
propagationDistance=100; %instead propagate all rays after they leave the drop by this amount


%set ray start points: 
xStart=linspace(-Rd, Rd, 2*Rd*rayDensity) %full width of droplet,
numRays=length(xStart);
yStart=(Rd+5)*ones(numRays); %5 units above droplet. 

%set ray start directions: 
rDir=repmat([0, -1],numRays ,1 );

%Ray amplitude is between 0 and 1, through away rays that have less than:
minimumAmplitude=0.1;

%% Run calculation
    
    % Make the drop object:
    theDrop=drop(dropLocation, Rd, Ri, no, ni, vr, dropAngle);

    % An Array to store the Rays in: 
    Rays=[];
    
    %coverslip data (planes that define top and bottom): 
    cTop=[ -Rd,0, 1];
    cBottom=[-Rd-cThickness, 0, 1];

    % Send Rays through scene
    for ii=1:numRays
        %create ray object: 
        aRay=ray([xStart(ii), yStart(ii)], rDir(ii, :), nm);
        
        %this is the line that actually does the raytracing calculation:
        aRayTraced=traceRay(aRay, theDrop, nm, minimumAmplitude);
        
        %propagate to screen and store
        for t=aRayTraced
            %coverslip: 
            if useCoverslip
                [loc, norm]=rayPlaneIntersection(t, cTop);
                if not(isnan(loc))
                    t=propagate(t, loc);
                    [t, tRef] = refract (t, norm, nc) ;
                end
                [loc, norm]=rayPlaneIntersection(t, cBottom);
                if not(isnan(loc))
                    t=propagate(t, loc);
                    [t, tRef] = refract (t, norm, na) ;
                end
            end
            
            
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
    
    %% plot the output: 
    figure
    
    drawDrop(theDrop) 
    %if you wish to modify colors of the drop, modify drawDrop() in Drop.m  
    
    % Draw the rays.  (Modify to draw only selected rays.)
    for Ray=Rays
        drawRay(Ray, [0, 0, 0.5]);
        hold on
    end
    
    if useCoverslip
    plot([-2*Rd, 2*Rd], [-Rd, -Rd], 'c')
    plot([-2*Rd, 2*Rd], [-Rd-cThickness, -Rd-cThickness], 'c')
    end
    
    axis image
    xlabel('x')
    ylabel('y')
    set(gcf, 'color', 'white')
