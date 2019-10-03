 function [ outRays ] = traceRayTriPhase( r, d, nm, minAmp, maxInteraction )
            %finds the next intersection points and refracts the ray.  
            tol=1e-14; %  <---how close two numbers should be to be considered equal
            if nargin<4
                minAmp=0.01; %Minimum intensity of ray to continue calculation (1% by default). 
            end
            if nargin <5
               maxInteraction=1000;  %Maximum times rays will reflect or refract within the drop. 
            end
            outRays=[];  
           
            %parameters of the drop
            Rd=d.radius;
            Ri1=d.innerRadius1;
            Ri2=d.innerRadius2;
            rot=d.rotAngle;
            
            %location of the innerDrop:
            if (length(d.location)==2) %2D
                p1Center=d.location+[d.d1*sind(rot), d.d1*cosd(rot)];
                p2Center=d.location+[d.d2*sind(rot), d.d2*cosd(rot)];
            end
            if (length(d.location)==3) %3D
                p1Center=d.location+[d.d1*sind(rot), 0 ,d.d1*cosd(rot)];
                p2Center=d.location+[d.d2*sind(rot), 0 ,d.d2*cosd(rot)];
            end
            
            %If the ray has already interacted with the drop too many
            %times, return empty.  Ray intersected too many times. 
            if(r.numIntersections>maxInteraction)
                %message='TOO MANY REFLECTIONS';
                return
            end
            
            %Ray on surface of drop (newly reflected Ray for example.)
            if(sum((r.location-d.location).^2) < Rd^2+tol && sum((r.location-d.location).^2) > Rd^2-tol && (r.location-d.location)*r.direction'> 0)
               if(abs(r.amplitude)>minAmp)
                   outRays=[r, outRays];
               end
               return %Ray is on surface of drop, facing outwards. 
            end
            
           
            
            
            %Is the ray outside the drop?
            if(sum((r.location-d.location).^2) > Rd^2+tol)
               
                %Will the ray hit the drop? if not exit
               [loc, norm]= raySphereIntersection( r, Rd, d.location);
               if (not(isnan(loc(1))))
                    r=propagate(r, loc);  %propagate the ray to that point
                    %Is this interface the inner or outer drop?
                    if (sum((loc-p1Center).^2) < Ri1^2+tol)
                        [r, rRef]=refract(r, norm, d.index1); 
                    elseif (sum((loc-p2Center).^2) < Ri2^2+tol) || ((sum((loc-p2Center).^2) > Ri2^2+tol)&& Ri2<0)
                        [r, rRef]=refract(r, norm, d.index2); 
                    else
                        [r, rRef]=refract(r, norm, d.index3); %outerDrop
    
                    end
                    %If the reflected Ray has high enough amplitude, trace
                    %it through the drop as well. 
                    if abs(rRef.amplitude)>=minAmp
                        rRef=traceRayTriPhase(rRef, d, nm, minAmp, maxInteraction);
                        outRays=[outRays, rRef];
                    end
                    
               else
                   outRays=r;
                   return %Doesn't hit the drop at all. 
               end
            end
            
            
            %Was enough of the ray transmitted, if not, exit. 
            if abs(r.amplitude)<minAmp
                return
            end
            

            %Do this up to four times: 
            for ii=1:4
                
            %Will the ray hit the interface 1 or interface 2, if both, which first?   
            [loc1, norm1]= raySphereIntersection( r, Ri1, p1Center);
            [loc2, norm2]= raySphereIntersection( r, abs(Ri2), p2Center);
            if isnan(loc1(1)) && isnan(loc2(1))
                break
            end
            if (not(isnan(loc1(1))) && sum((loc1-d.location).^2)<Rd^2)&& (not(isnan(loc2(1))) && sum((loc2-d.location).^2)<Rd^2)
                if sum((loc1-r.location).^2)<sum((loc2-r.location).^2)
                    loc2(:)=nan;
                else
                    loc1(:)=nan;
                end
            end
            
            if(not(isnan(loc1(1))) && sum((loc1-d.location).^2)<Rd^2)  %Intersection found and is inside full drop.

                r=propagate(r, loc1);  %propagate the ray to that point
     
                if (r.direction*(r.location-p1Center)' >=0 ) %Ray pointing inwards or outwards?
                    [r, rRef]=refract(r, norm1, d.index2);  %Ray was facing outwards
                else
                    [r, rRef]=refract(r, norm1, d.index1); %Ray was facing inwards. 
                end
                    if abs(rRef.amplitude)>=minAmp
                        rRef=traceRayTriPhase(rRef, d, nm, minAmp, maxInteraction);
                        outRays=[outRays, rRef];
                    end
                
            end
            if(not(isnan(loc2(1))) && sum((loc2-d.location).^2)<Rd^2)  %Intersection found and is inside full drop.
                r=propagate(r, loc2);  %propagate the ray to that point
                if ((r.direction*(r.location-p2Center)' >=0 && Ri2>0) || (r.direction*(r.location-p2Center)' <=0 && Ri2<0) ) %Ray pointing inwards or outwards?
                    [r, rRef]=refract(r, norm2, d.index3);

                else
                    [r, rRef]=refract(r, norm2, d.index2);

                end
                    if abs(rRef.amplitude)>=minAmp
                        rRef=traceRayTriPhase(rRef, d, nm, minAmp, maxInteraction);
                        outRays=[outRays, rRef];
                    end
                
            end
                        %Was enough of the ray transmitted, if not, exit.
            if abs(r.amplitude)<minAmp
                return
            end
            end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
%             %Still inside phase 1?  Find next outward interface. 
%             if (r.direction*(r.location-p1Center)' <= 0)  && sum((r.location-p1Center).^2<=Ri1^2+tol) 
%                     [loc, norm]= raySphereIntersection( r, Ri1, p1Center);
%                     if sum((loc-d.location).^2)<Rd^2+tol && not(isnan(loc(1)))
%                         r=propagate(r, loc);
%                         [r, rRef]=refract(r, -norm, d.index2);
%                         if abs(rRef.amplitude)>=minAmp
%                             rRef=traceRayTriPhase(rRef, d, nm, minAmp, maxInteraction);
%                             outRays=[outRays, rRef];
%                         end
%                     end
%             end
%             
%             
%             %Still inside phase 3?  Find next outward interface. 
%              inPhase3=(r.direction*(r.location-p2Center)' <= 0  && sum((r.location-p2Center).^2<=Ri2^2+tol) && Ri2<0)...
%                 ||not((r.direction*(r.location-p2Center)' <= 0  && sum((r.location-p2Center).^2<=Ri2^2+tol) && Ri2>0));
%             if false && inPhase3
%                     [loc, norm]= raySphereIntersection( r, abs(Ri2), p2Center);
%                     if sum((loc-d.location).^2)<Rd^2+tol && not(isnan(loc(1)))
%                         r=propagate(r, loc);
%                         [r, rRef]=refract(r, -norm, d.index2);
%                         if abs(rRef.amplitude)>=minAmp
%                             rRef=traceRayTriPhase(rRef, d, nm, minAmp, maxInteraction);
%                             outRays=[outRays, rRef];
%                         end
%                     end
%             end  
%             
%             %Inside phase 2?  Find next outward interface. 
%             inPhase1=(r.direction*(r.location-p1Center)' <= 0  && sum((r.location-p1Center).^2<=Ri1^2+tol))
%             %inPhase3=((r.direction*(r.location-p2Center)' <= 0  && sum((r.location-p2Center).^2<=Ri2^2+tol) && Ri2<0)||...
%             
%             inPhase3=not((r.direction*(r.location-p2Center)' <= 0  && sum((r.location-p2Center).^2<=Ri2^2+tol) && Ri2>0))
%             if not(inPhase1 || inPhase3)
%                 test='here'
%                     [loc, norm]= raySphereIntersection( r, abs(Ri2), p2Center);
%                     if sum((loc-d.location).^2)<Rd^2+tol && not(isnan(loc(1)))
%                         r=propagate(r, loc);
%                         [r, rRef]=refract(r, -norm, d.index3);
%                         if abs(rRef.amplitude)>=minAmp
%                             rRef=traceRayTriPhase(rRef, d, nm, minAmp, maxInteraction);
%                             outRays=[outRays, rRef];
%                         end
%                     else    
%                         [loc, norm]= raySphereIntersection( r, Ri1, p1Center);
%                     if sum((loc-d.location).^2)<Rd^2+tol && not(isnan(loc(1)))
%                         r=propagate(r, loc);
%                         [r, rRef]=refract(r, -norm, d.index1);
%                         if abs(rRef.amplitude)>=minAmp
%                             rRef=traceRayTriPhase(rRef, d, nm, minAmp, maxInteraction);
%                             outRays=[outRays, rRef];
%                         end
%                         
%                     end
%                     end
%             end
          

   
            %%
            
            %Was enough of the ray transmitted, if not, exit.
            if abs(r.amplitude)<minAmp
                return
            end
            
            %Ray hits outside interface of drop:
            [loc, norm]= raySphereIntersection( r, Rd, d.location);  %Find final intersection
            if(not(isnan(loc(1)))) %Check that it exists. 
                r=propagate(r, loc);  %propagate the ray to that point
                [r, rRef]=refract(r, -norm, nm); %And refract at that interface. 
                if abs(rRef.amplitude)>=minAmp
                    rRef=traceRayTriPhase(rRef, d, nm, minAmp, maxInteraction);
                    
                    outRays=[outRays, rRef];
                end
                
            end
            if abs(r.amplitude)>=minAmp
            outRays=[outRays, r]; %Add the ray to the array. 
            end
       end


