%There are three primary components of the raytracer: 
%rays, drops, and the traceRay function. 

%The raytracer is untiless, but easy enough to assume all units are in
%micrometers.

% Make sure DropletTracer folder is in your path. 
addpath('C:\Users\Sara\Dropbox (MIT)\DropletRayTracer-Full\DropletTracer')

%% Droplet parameters:
Rd=50
Ri1=50%00000
Ri2=50%00000
nm=1.33
ind1=1.38
ind2=1.41
ind3=1.27
v1=0.3
v2=0.4
minimumAmplitude=0.5


%% Raytracing Parameters:
%set number of rays: 
rayDensity=0.1 %rays per unit length

%include coverslip (or any glass substrate)
useCoverslip=false;
if useCoverslip
    cThickness=105%thickness of coverslip
    nc=1.52 %refractive index of coverslip
    na=1% refractive index outside of coverslip (air probably)
end
%if you use a thicker substrate, check the location of the screen so that it is after the bottom edge of the coverslip
screenDistance=Rd+500;
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
    theDrop=TriPhaseDrop([0, 0], Rd, Ri1, Ri2, ind1, ind2, ind3, v1, v2)

    % An Array to store the Rays in: 
    Rays=[];
    
    %coverslip data (planes that define top and bottom): 
    if useCoverslip
        cTop=[ -Rd,0, 1];
        cBottom=[-Rd-cThickness, 0, 1];
    end

    % Send Rays through scene
    for ii=1:numRays
        %create ray object: 
        aRay=ray([xStart(ii), yStart(ii)], rDir(ii, :), nm);
        
        %this is the line that actually does the raytracing calculation:
        aRayTraced=traceRayTriPhase(aRay, theDrop, nm, minimumAmplitude);
        
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
    % Draw the rays.  (Modify to draw only selected rays.)
    for Ray=Rays
        drawRay(Ray);
        hold on
    end
       
    drawDrop(theDrop) 
    %if you wish to modify colors of the drop, modify drawDrop() in Drop.m  
    
    if useCoverslip
    plot([-2*Rd, 2*Rd], [-Rd, -Rd], 'c')
    plot([-2*Rd, 2*Rd], [-Rd-cThickness, -Rd-cThickness], 'c')
    end
    
    axis image
    xlabel('x')
    ylabel('y')
    set(gcf, 'color', 'white')
