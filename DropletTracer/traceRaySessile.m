function [ outRays ] = traceRaySessile( r, Rd, CA, n, nm, minAmp, maxInteraction )
            %finds the next intersection points and refracts the ray.  
            % r: ray object
            %Rd: radius of the sessile drop
            %CA: contact angle of the drop (in degrees)
            % n: refractive index of drop
            %nm: refractive index of surrounding medium
            %minAmp: minimum Amplitude to keep ray
            %maxInteraction: maximum number of reflection/refraction events
            %before throwing away ray (in case of trapped rays.)
            %center of curvature is located at [0,0, 0]
            
            
            tol=1e-8; %  <---how close two numbers should be to be considered equal
            
            if nargin<5
                nm=1 %Refractive index of surrounding medium
            end
            
            if nargin<6
                minAmp=0.01; %Minimum intensity of ray to continue calculation (1% by default). 
            end
            if nargin <7
               maxInteraction=1000;  %Maximum times rays will reflect or refract within the drop. 
            end
            outRays=[];  
           
            d=-cosd(CA)*Rd %distance from center of curvature to interface
            
            %plane defining substrate: 
            if length(r.location)==2
                dim=2;
                substrate=[d, 0, 1];
            elseif length(r.location)==3
                dim=3;
                substrate=[d, 0, 0, 1];
            end
            
            
            
            %If the ray has already interacted with the drop too many
            %times, return empty.  Ray intersected too many times. 
            if(r.numIntersections>maxInteraction)
                %message='TOO MANY REFLECTIONS';
                return
            end
            
            %Ray on surface of drop (newly reflected Ray for example.)
            if(sum((r.location).^2) < Rd^2+tol && sum(r.location.^2) > Rd^2-tol && r.location(end)<d && (r.location)*r.direction'> 0)
               if(abs(r.amplitude)>minAmp)
                   outRays=[r, outRays];
               end
               return %Ray is on surface of drop, facing outwards. 
            end
            
            %ray above substrate: 
            if (r.location(end) > d)
                message='section 1'
                %Find intersection with substrate
                [loc, norm]= rayPlaneIntersection(r, substrate);
                                    
                if (not(isnan(loc(1))))
                    r=propogate(r, loc);
                    if(sum((r.location).^2) < Rd^2+tol)
                    [r, rRef]=refract(r, norm, n);
                     ramp=r.amplitude
                    if abs(rRef.amplitude)>=minAmp
                        outRays=[outRays, rRef];
                    end
                    end
                end
                ramp=r.amplitude
            end
            
            %Is the ray outside the drop, but at or below substrate?
            if(sum((r.location).^2) > Rd^2+tol && r.location(end)< d+tol)
               %Will the ray hit the drop? if not exit
               message='section 2'
               currentLoc=r.location
               [loc, norm]= raySphereIntersection( r, Rd, 0*(1:dim))
               if (not(isnan(loc(1))))
                    r=propogate(r, loc);  %Propogate the ray to that point
                    [r, rRef]=refract(r, norm, nm); %outerDrop
                    %If the reflected Ray has high enough amplitude, trace
                    %it through the drop as well. 
                    if abs(rRef.amplitude)>=minAmp
                        rRef=traceRaySessile(rRef, Rd, CA, n, nm, minAmp, maxInteraction);
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
            
           
            %Ray hits outside interface of drop:
            [loc, norm]= raySphereIntersection( r, Rd, 0*(1:dim));  %Find final intersection
            if not(isnan(loc(1))) %Check that it exists. It should.
                if loc(end)< d
                    currentLoc=loc
                    currentD=d
                    ramp=r.amplitude
                    r=propogate(r, loc);  %Propogate the ray to that point
                    currentNorm=norm
                    [r, rRef]=refract(r, norm, nm); %And refract at that interface. 
                    ramp=r.amplitude
                elseif (loc(end)>d)
                    [loc, norm]= rayPlaneIntersection(r, substrate)
                    r=propogate(r, loc);  %Propogate the ray to that point
                    [r, rRef]=refract(r, norm, nm); %And refract at that interface. 
                end

                            
            if abs(rRef.amplitude)>=minAmp
                rRef=traceRaySessile(rRef, Rd, CA, n, nm, minAmp, maxInteraction);
                outRays=[outRays, rRef];
            end
            end
            if abs(r.amplitude)>=minAmp
            outRays=[outRays, r]; %Add the ray to the array. 
            end
       end


