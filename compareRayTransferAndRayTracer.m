clear all
close all


%There are three primary components of the raytracer: 
%rays, drops, and the traceRay function. 

%The raytracer is untiless, but easy enough to assume all units are in
%micrometers.

% Make sure DropletTracer folder is in your path. 
addpath('DropletTracer')

%% Droplet parameters:
vr=1; %volume ratio.  V_internal/v_external
Rd=1; %droplet radius
dropLocation=[0,0];
dropAngle=0;  % rotation angle (degrees)
Ri=1; %Internal radius of curvature. use positve numbers

% Indexes of refraction: 
nm=1.33; % Surrounding Medium
no=1.27; % Outer Phase
ni=1.39;% Inner Phase

offAxisDists=[0.7, 0.5, 0.3, 0.1]

%% Raytracing Parameters:
%set number of rays: 
rayDensity=5 %rays per unit length

%include coverslip (or any glass substrate)
useCoverslip=false;
if useCoverslip
    cThickness=105%thickness of coverslip
    nc=1.52 %refractive index of coverslip
    na=1% refractive index outside of coverslip (air probably)
end
%if you use a thicker substrate, check the location of the screen so that it is after the bottom edge of the coverslip
screenDistance=Rd+10;
%plane to propagate rays to after interaction with the droplet
%if rays don't hit this plane they are thrown away:
%equation of plane is p(1)=p(2)x+p(3)y. 
screen=[ -screenDistance,0, 1]; % Plane that will be the screen.  
screenSize=200
%you can choose to not use a screen: 
useScreen=true; 
propagationDistance=100; %instead propagate all rays after they leave the drop by this amount


%set ray start points: 
xStart=[-offAxisDists, 0] %full width of droplet,
numRays=length(xStart);
yStart=(Rd+0.5)*ones(numRays); %5 units above droplet. 

%set ray start directions: 
rDir=repmat([0, -1],numRays ,1 );

%Ray amplitude is between 0 and 1, through away rays that have less than:
minimumAmplitude=0.1;

%% Display Rays for one geometry: 
    
    % Make the drop object:
    theDrop=drop(dropLocation, Rd, Ri, no, ni, vr, dropAngle);

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
            [loc, norm]=rayPlaneIntersection(t, [0, 1, 0]);
            %if (abs(loc(1))< screenSize)
                t=propagate(t, loc);
                Rays=[Rays, t];
            %end
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
    xlim([-1.5, 1.5])
    axis image
    xlabel('x')
    ylabel('y')
    set(gcf, 'color', 'white')
    
    
    
    %%
    
        aRay=ray([xStart(1), yStart(1)], rDir(1, :), nm);
        
        %this is the line that actually does the raytracing calculation:
        aRayTraced=traceRay(aRay, theDrop, nm, minimumAmplitude)
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
        end
        if length(t)>1
            t=t(1)
        end
        %where does t hit the optical axis: 
        [loc, norm]=rayPlaneIntersection(t, [0, 1, 0]);
            
%%  Plot Comparing Ray Tracer and Ray Transfer: 
Ris=0.8:0.01:3
nH=ni
nF=no
for jj=1:length(Ris)
    for kk=1:length(xStart)
    Ri=Ris(jj);
    theDrop=drop(dropLocation, Rd, Ri, no, ni, vr, dropAngle);
    aRay=ray([xStart(kk), yStart(1)], rDir(1, :), nm);
        
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
        end
        if length(t)>1
            t=t(1)
        end
        %where does t hit the optical axis: 
        [loc, norm]=rayPlaneIntersection(t, [0, 1, 0]);
        BFL_Tr(jj, kk)=-1-loc(2);
    end 
        
        
        
        %Ray Transfer Matrices
    d=theDrop.offCenterDist; 
    l1 = Rd - d + Ri;
    l2 = Rd - Ri + d;

    intoDrop = [[1,  0]; [(nm - nH)/(Rd*nH), nm/nH]];
    throughFirstPhase = [[1, l1]; [0, 1]];
    intoInterface = [[1,  0]; [(nH - nF)/(-Ri*nF), nH/nF]];
    throughSecondPhase = [[1, l2]; [0, 1]];
    outOfDrop = [[1,  0]; [(nF - nm)/(-1*Rd*nm), nF/nm]];

    dropMatrix=outOfDrop*throughSecondPhase*intoInterface*throughFirstPhase*intoDrop;
    
    BFL_m(jj)=-dropMatrix(1,1)/dropMatrix(2,1);
end
%%
figure
plot(Ris, BFL_m, 'r', 'linewidth', 2)
hold on
plot(Ris, BFL_Tr(:, 4), 'linewidth', 1, 'color', [0, 0, 1])
plot(Ris, BFL_Tr(:, 3), 'linewidth', 1, 'color', [0, 0.5, 1])
plot(Ris, BFL_Tr(:, 2), 'linewidth', 1, 'color', [0, 0.7, 1])
plot(Ris, BFL_Tr(:, 1), 'linewidth', 1, 'color', [0, 1, 1])

xlabel('R_i/R_d')
ylabel('BFL/R_d')
set(gcf, 'color', 'w')
legend('Ray Transfer Matrix', 'Ray Tracing x_i=0.1','Ray Tracing x_i=0.3','Ray Tracing x_i=0.5','Ray Tracing x_i=0.7' , 'location', 'best')










