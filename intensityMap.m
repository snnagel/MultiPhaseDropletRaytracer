clear all
close all

resolution=20;
dropRadius=1;
rayDensity=resolution*5; %start with at least 5 rays per pixel, spread out evenly
rayStartHeight=1;  %max height of starting ray (spread out evenly from +-0)
dropCenter=[2, 0];


normalization=(resolution/rayDensity)^2;  %intensity per pixel before the drop. 

%Indexes of refraction: 
nm=1.33; %surrounding Medium
ni=1.387;%innerPhase
no=1.27; %outerPhase
vr=0.45

screenHeight=2; %The height in the x2 direction that will be mapped. 
distance=20;

Ri=(vr/(1+vr))^(1/3);

startPlane=[dropCenter(1)+dropRadius, 1, 0]; %After exiting the drop, propogate rays to this point before measuring intensity. 

theDrop=drop(dropCenter, dropRadius, Ri, no, ni, vr, 180);

    %An Array to store the Rays in: 
    Rays=[];

    %Send Rays through scene
    for x=[-rayDensity*rayStartHeight:rayDensity*rayStartHeight]
        c=ray([0, x/rayDensity], [1, 0], nm);
        cTraced=traceRay(c, theDrop, nm)  %Trace Ray through the drop
        for t=cTraced
            [loc, norm]=rayPlaneIntersection(t, startPlane); %Propogate to the start plane
            t=propogate(t, loc)
            Rays=[Rays, t];
        end
    end
    
    %Draw the rays:
    figure
    for ray=Rays
        drawRay(ray)
        hold on
    end
    drawnow
    drawDrop(theDrop)
    axis equal
    
    
    %%
    numSlices=resolution*distance;
    Intensity=zeros(screenHeight*2*resolution, numSlices);
    for k=[1:numSlices]
        I=IntensityPlane(Rays, startPlane(1)+k/resolution, resolution, screenHeight);
        Intensity(:, k)=I/normalization;
    end
    %%
    figure

    
    
 
    for ray=Rays
        ray=propogate(ray, distance);
        drawRay(ray)
        hold on
    end
    drawDrop(theDrop)
    axis equal
    
    
    z=[0:distance]+startPlane(1);
    x=[-screenHeight:screenHeight];
    imagesc(z, x, Intensity)
    axis image
    colorbar
    
    axis([0, max(z), -screenHeight screenHeight])
    
   

