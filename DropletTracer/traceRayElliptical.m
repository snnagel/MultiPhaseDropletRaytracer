function [ outRays ] = traceRayElliptical( r, d, nm, minAmp, maxInteraction )
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
            a=d.a;
            b=d.b;
            
            rot=d.rotAngle;
            
            %location of the innerDrop:
            if (length(d.location)==2) %2D
                innerDropCenter=d.location+[d.offCenterDist*sind(rot), d.offCenterDist*cosd(rot)];
            end
            if (length(d.location)==3) %3D
                innerDropCenter=d.location+[d.offCenterDist*sind(rot), 0 ,d.offCenterDist*cosd(rot)];
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
                    r=propogate(r, loc);  %Propogate the ray to that point
                    %Is this interface the inner or outer drop?
                    if (sum((loc(1:end-1)-innerDropCenter(1:end-1)).^2/b^2+(loc(end)-innerDropCenter(end))^2/b^2 < 1+tol))
                        [r, rRef]=refract(r, norm, d.innerIndex); %innerDrop
                    else
                        [r, rRef]=refract(r, norm, d.outerIndex); %outerDrop
                    end
                    %If the reflected Ray has high enough amplitude, trace
                    %it through the drop as well. 
                    if abs(rRef.amplitude)>=minAmp
                        rRef=traceRayElliptical(rRef, d, nm, minAmp, maxInteraction);
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
            

            
            %Will the ray hit the inside drop?              
            [loc, norm]= rayEllipseIntersection( r, a, b, innerDropCenter); %Try to find intersection with innerDrop
            if(not(isnan(loc(1))) && sum((loc-d.location).^2)<Rd^2)  %Intersection found and is inside full drop.

                r=propogate(r, loc);  %Propogate the ray to that point
                if (r.direction*(norm)'>=0) %Ray pointing inwards or outwards?
                    [r, rRef]=refract(r, norm, d.outerIndex);  %Ray was facing outwards
                else
                    [r, rRef]=refract(r, norm, d.innerIndex); %Ray was facing inwards. 
                    
                end
                    if abs(rRef.amplitude)>=minAmp
                        rRef=traceRayElliptical(rRef, d, nm, minAmp, maxInteraction);
                        outRays=[outRays, rRef];
                    end
                
            end
            
            
            
            %Was enough of the ray transmitted, if not, exit.
            if abs(r.amplitude)<minAmp
                return
            end
            
                %Still inside innerDrop?  Find next outward interface. 
                if length(r.location)==2  %2D case: 
                    insideInner=(r.location(1)-innerDropCenter(1)).^2/b^2+(r.location(2)-innerDropCenter(2))^2/a^2 < 1+tol;
                    facingInwards=r.direction*[r.location(1)/b^2, r.location(2)/a^2]';
                elseif length(r.location)==3 % %3D case
                    insideInner= sum((r.location(1:2)-innerDropCenter(1:2)).^2)/b^2+(r.location(3)-innerDropCenter(3))^2/a^2 < 1+tol;
                    facingInwards=r.direction*[r.location(1)/b^2, r.location(2)/b^2, r.location(3)/a^2]'<0;
                end
                %make sure that only true interface is hit: 
                belowContactPoint=(r.location(end)<(d.interfaceZ-d.location(end)));
                
                
            if (facingInwards  && insideInner && belowContactPoint) 
                    [loc, norm]= rayEllipseIntersection( r, a, b, innerDropCenter);
                    if sum((loc-d.location).^2)<Rd^2+tol && not(isnan(loc(1)))
                        r=propogate(r, loc);
                        [r, rRef]=refract(r, -norm, d.outerIndex);
                        if abs(rRef.amplitude)>=minAmp
                            rRef=traceRayElliptical(rRef, d, nm, minAmp, maxInteraction);
                            outRays=[outRays, rRef];
                        end
                    end        
            end
            
            
            %Was enough of the ray transmitted, if not, exit.
            if abs(r.amplitude)<minAmp
                return
            end
            
            %Ray hits outside interface of drop:
            [loc, norm]= raySphereIntersection( r, Rd, d.location);  %Find final intersection
            if(not(isnan(loc(1)))) %Check that it exists. 
                r=propogate(r, loc);  %Propogate the ray to that point
                [r, rRef]=refract(r, -norm, nm); %And refract at that interface. 
                
                if abs(rRef.amplitude)>=minAmp
                    rRef=traceRayElliptical(rRef, d, nm, minAmp, maxInteraction);
                    outRays=[outRays, rRef];
                end
                
            end
            if abs(r.amplitude)>=minAmp
            outRays=[outRays, r]; %Add the ray to the array. 
            end
       end


