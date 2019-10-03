function IntMap = IntensityPlaneVert( rayArray, planeLoc, resolution, screenSize )
%Takes a bundle of rays (rayArray) and adds up Intensity INCOHERENTLY through a plane
%normal to x1 axis.  with location planeLoc
%%resolution is the number of pixels per unit size. 

numPixels=resolution*screenSize*2;  %total number of pixels across. 

%First initialize screen and plane

%2D:
if length(rayArray(1).location)==2
    plane=[planeLoc, 0, 1];
    screen=zeros(numPixels, 1);
    dim=2;
end

%3D
if length(rayArray(1).location)==3
    plane=[planeLoc, 0, 0, 1];
    screen=zeros(numPixels);
    dim=3;
end


%look at each ray through the screen:
for ray=rayArray
    [loc, norm]=rayPlaneIntersection(ray,plane); %see ray.m

        pixNum=round((loc(1:end-1)+screenSize)*resolution);   %Round to nearest whole pixel number
        if pixNum > 0 && pixNum <= numPixels %Hits the inside the screen 
        
            % add the amplitude of the ray to that pixel. 
        if dim==2 %2D case
            screen(pixNum)=screen(pixNum)+abs(ray.amplitude)^2;
        end
        if dim==3 %3D case
            screen(pixNum(1), pixNum(2))=screen(pixNum(1), pixNum(2))+abs(ray.amplitude)^2;
        end 
    end
end 
IntMap=screen;
    
       
end

