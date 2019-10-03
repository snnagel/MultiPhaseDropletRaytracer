classdef ray
    properties 
        location	%Current location of the ray (z, x)
		direction	%Current direction of the ray (Unit vector (z, x)
		trajectory	%Array of past locations of the ray
        numIntersections %number of intersections that ray has had
        currentIndex
        amplitude  %Intensity is abs(amplitude)^2. note that this is actually the |amp| phase is included in polarization
        amplitudeHistory  %Array of past amplitudes of the ray.  
        polarization %Polarization vector. Points in direction of electric field, and includes amplitude and phase information
        TIR=false %tag to say whether at any time total internal reflection occured. 
        OPL %optical path length- updated when propagate is called.
    end
    
    methods
		function r=ray(loc, dir, n, amp, pol) %Constructor  
            %The location and direction are required, index and intensity
            %are optional (default 1).  If intensity is specified, so must the index. 
			r.location=loc;
            r.direction=dir/sqrt(sum(dir.^2));  %Direction must be normalized. unit vector. 
            r.trajectory=loc;  %Add virst location to trajectory list. 
            r.numIntersections=0; 
            
            if nargin>2
                r.currentIndex=n;
            else
                r.currentIndex=1;
            end
            if nargin >3
                r.amplitude=amp;
            else 
                r.amplitude=1;
            end
            r.amplitudeHistory=r.amplitude;
            
            if nargin < 5
                %polarization and direction must be perpendicular: 
                pol=[0, 0, 1];
            end
            if length(dir)==3 
                if sum(dir==pol)==3 || sum(dir==-pol)==3
                pol=rand(1, 3); %Assigns random polarization. 
                end
            end
            if length(dir)==2 
                if sum([dir, 0]==pol)==3 || sum([dir, 0]==-pol)==3
                pol=rand(1, 3); %Assigns random polarization. 
                end
            end
            if length(dir)==3
                pol=pol-dir*pol'*dir;
            else
                pol=pol-[dir, 0]*pol'*[dir, 0];
            end
                r.polarization=r.amplitude*pol/sqrt(pol*pol');
            
            r.OPL=0;
                       
        end
        
        %Function to set direction. 
        function r = set.direction(r,dir)
            r.direction=dir/sqrt(dir*dir');  %direction must always have unit length. 
        end
        
        function r=set.amplitude(r, amp)
            r.amplitude=amp;
            r.amplitudeHistory=[r.amplitudeHistory, amp];
        end
            
        %propagate changes the current location of the ray, and updates the
        %trajectory.  Later might want to add a phase change associated
        %with this. Currently doesn't check that the ray CAN go to that
        %location..
        %If only a scala is specified the ray travels that many units in
        %the direction it currently points. 
		function rOut = propagate(r, l)		
			%update location and trajectory
            if size(l)==1 %1D value for the propogation given, go l units in current ray direction. 
                r.location=r.location+l.*r.direction;
                r.OPL=r.OPL+l*r.currentIndex;
            else
                tl=sqrt(sum((l-r.location).^2));
                r.location=l;
                r.OPL=r.OPL+tl*r.currentIndex;
            end
			r.trajectory=[r.trajectory; r.location];
            rOut = r;
        end
        
        
        
        %Refract a ray through interface with normal norm and index n2. 
		function [rayTransmitted, rayReflected] = refract (r, norm, n2) 
            
            %2D or 3D?
            dim=length(r.location);
            
            n1=r.currentIndex;
            rn=n1/n2; %Index ratio
            
            cosi=-norm*r.direction'; %cosine of Incoming angle
            
            r.numIntersections=r.numIntersections+1;
            
            
            if cosi<0
                norm=-norm;
                cosi=-cosi;
            end

            l=r.direction;

            %snells law:
            s=1-rn^2*(1-cosi^2);
            v=rn*l+(rn*cosi-sqrt(s))*norm; %Vector Form of snell's law. 
            
            %Fresnel Equations:
            %cost=v*norm'/sqrt(v*v'); %cosine of Transmitted angle
            cost=sqrt(1-(n1/n2)^2*(1-cosi^2));
            if s>0
            
            rp=(n2*cosi-n1*cost)./(n2*cosi+n1*cost);
            rs=(n1*cosi-n2*cost)/(n1*cosi+n2*cost);
            ts=(2*n1*cosi)/(n1*cosi+n2*cost);
            tp=(2*n1*cosi)/(n1*cost+n2*cosi);
            
            else
                ts=0;
                tp=0;
                n21=n2/n1;
                rs=(cosi-1i*sqrt(1-cosi^2-n21^2))/(cosi+1i*sqrt(1-cosi^2-n21^2));
                rp=(n21^2*cosi-1i*sqrt(1-cosi^2-n21^2))/(n21^2*cosi+1i*sqrt(1-cosi^2-n21^2));
            end
            
            %Need to consider Polarization:
            %initial polarization:
            pol=r.polarization;
            
            %normal to plane of incidence:
           if dim==2
               PlaneOfIncidence=[0, 0 , 1];
           end
           if dim==3
               PlaneOfIncidence=cross(norm, r.direction);
               %normalize it:
               PlaneOfIncidence=PlaneOfIncidence/sqrt(PlaneOfIncidence*PlaneOfIncidence');
           end
           
           %normal incidence, treat as TE, same thing.
           if(cosi==1)
                if dim==2
                   sPol=pol;
                   pPol=0*cross(sPol, [norm, 0]);
                   %pPolAmp=0;
                   PlaneOfIncidence=sPol;
                end
                if dim==3
                   sPol=pol;
                   pPol=0*cross(sPol, norm);
                   %pPolAmp=0;
                   PlaneOfIncidence=sPol;
               end
           else
           PlaneOfIncidence=PlaneOfIncidence/sqrt(PlaneOfIncidence*PlaneOfIncidence');
           
           %TE polarization component:
           sPol=(pol*PlaneOfIncidence')*PlaneOfIncidence;
           
           %TM component:
           pPol=pol-sPol;

           end
           
           refDir=r.direction+2*cosi*norm;
           refDir=refDir/sqrt(refDir*refDir');
           
           if dim==2
               dirR3=[refDir, 0];
               dirT3=[v, 0];
           else
               dirR3=refDir;
               dirT3=v;
           end
           
        %_______________________________________________________   
           %rotate the TM polarization: 
           k=PlaneOfIncidence; %unit vector to rotate around
           thetar=acos(cosi);
   
           pPolRotr=pPol*cos(2*thetar)+cross(k, pPol)*sin(2*thetar)+k*(k*pPol')*(1-cos(2*thetar));
           if not(abs(pPolRotr*dirR3')<1e-14)
               %'Rotated Wrong Way'
               pPolRotr=pPol*cos(-2*thetar)+cross(k, pPol)*sin(-2*thetar)+k*(k*pPol')*(1-cos(-2*thetar));
           
           end
           
           thetat=acos(cosi)-acos(cost);
           pPolRott=pPol*cos(thetat)+cross(k, pPol)*sin(thetat)+k*(k*pPol')*(1-cos(thetat));
           if s>0
           if not(abs(pPolRott*dirT3')<1e-14)
               %'Rotated Wrong Way'
               pPolRott=pPol*cos(-thetat)+cross(k, pPol)*sin(-thetat)+k*(k*pPol')*(1-cos(-thetat));
           
           end
           polTrans=tp*pPolRott+ts*sPol;
           else
               polTrans=[0, 0, 0];
           end
           
           polRef=rp*pPolRotr+rs*sPol;

          %_______________________________________________________  
                    
          

            rayT=r;
            rayT.currentIndex=n2; 
            rayT.direction=v/sqrt(v*v');
            rayT.polarization=polTrans;
            rayT.amplitude=sqrt(polTrans*polTrans');

            rayR=r;
            rayR.direction=refDir;
            rayR.polarization=polRef;
            rayR.amplitude=sqrt(polRef*polRef');
            if s<=0
                rayR.TIR=true;
                rayT.amplitude=0;
            end
            
            %The "History" of the ray will be held in both of the rays,
            %Beware of duplication if summing over rays wthin structure...

              rayTransmitted=rayT;
              rayReflected=rayR;
            
        end
        
        
        
        %reflect a ray from surface with normal, norm. reflectance is the
        %ratio of amplitude reflected. (the rest is either tranmissited or
        %absorbed. defaults to 1. 
        function rRef= reflect( r, norm, reflectance)
            cosi=norm*r.direction'; %Incoming angle
            v=r.direction-2*cosi*norm;
            rRef=r;
            rRef.direction=v;
            %if user specified a reflectance, reduce intensity by that
            %fraction. 
            if(nargin ==3)
                rRef.amplitude=rRef.amplitude*reflectance;
            end
        end
        
        
        
        function drawRay(r, rcolor)
            if nargin<2 
                rcolor='k'
            end
               
            %Draw in 2D
            if(length(r.location)==2)
                plot(r.trajectory(:,1), r.trajectory(:, 2), 'color', rcolor, 'Linewidth', 1.5)              
            end
            
            %Or in 3D:
            if(length(r.location)==3)
                plot3(r.trajectory(:,1), r.trajectory(:, 2), r.trajectory(:, 3), 'color', rcolor);
            end

        end
        
        
        
        function [ location, normal ] = raySphereIntersection( r, Rad, loc )
           %returns the loc and normal of intersection between a sphere of radius Rad with center at location loc and a ray r
            %ray position 
            tol=1e-8;
           
            %distance center of sphere from ray start point:
            dc=r.location-loc;
            
            %Components of equation to give intersection distances. 
            b=2*sum(dc.*r.direction);
            c=sum(dc.^2)-Rad^2;
            delta=b^2-4*c;
            
            location=nan;  %If location is not assigned a value, nan is returned. 
            if delta > tol
                u1 = (-b -sqrt(delta)) / 2;
                u2 = (-b +sqrt(delta)) / 2;
                
                %want the nearest intersection, but in the direction the
                %ray is pointing.  If the nearest intersection is within
                %the tolerance, assume the ray already saw it.  
                if ((u1 > u2 || u1 < tol) && u2 > tol  )  
                    location= r.location+u2*r.direction;
                end
                if ((u1 < u2 || u2 < tol )&& u1 > tol )
                   location= r.location+u1*r.direction;
                end 
            else
                if abs(delta) < tol
                % delta around zero: find unique root.
                u = -b/2;    
                location = r.location + u*r.direction;
                end
            end
            normal=(location-loc)/Rad; %unit normal at intersection point.  
        end
       
        %Calculates intersection point with a plane.  The plane is defined
        %by p:
        %p(2)*x1+p(3)*x2+p(4)*x3= p(1) is the equation of the plane. 
        function [location, normal]= rayPlaneIntersection(r, p)
            %normal vector to plane:
            %Not necessarily normalized yet, which is ok.  otherwise would
            %need to change p(1). 
            normal=p(2:end);
            
            if abs(normal*r.direction') > 10^(-14)  %make sure ray isn't parallel to screen
                %Distance along direction vector to intersection:
                u=(p(1)-normal*r.location')/(normal*r.direction');

                %if the intersection is backwards on ray trajectory ignore it. 
                if u<0
                    u=nan;
                end
            else 
                u=nan;
            end 
            %Return location of intersection.  
            location=r.location+u*r.direction;
            normal=normal/(normal*normal');  %Normalize now.  
        end
        
        
        
        
        
        
        
        function [ location, normal ] = rayEllipseIntersection( r, a, b, loc )
            %returns the loc and normal of intersection between an elliptoid with axes a and b centered at location loc with ray r
            %Orientation is the direction of the symmetry axis of the
            %elliptoid (a-direction)
            %ray position 
            tol=1e-8;
           
            %distance center of ellipse from ray start point:

            dc=r.location-loc;
            
            %Seperate 2D and 3D cases: 
            if length(dc) == 2
                xc=dc(1);
                zc=dc(2);
                xd=r.direction(1);
                zd=r.direction(2);
                
                %Parts of quadratic formula: 
                A=xd^2/b^2+zd^2/a^2;
                B=2*xc*xd/b^2+2*zc*zd/a^2;
                C=xc^2/b^2+zc^2/a^2-1;
  
            end
            if length(dc) == 3
                xc=dc(1);
                yc=dc(2);
                zc=dc(3);
                xd=r.direction(1);
                yd=r.direction(2);
                zd=r.direction(3);
                
                
                A=(xd^2+yd^2)/b^2+zd^2/a^2;
                B=2*(xc*xd+yc*yd)/b^2+2*zc*zd/a^2;
                C=(xc^2+yc^2)/b^2+zc^2/a^2-1;
                
            end

            delta=B^2-4*A*C;
            
            location=nan;  %If location is not assigned a value, nan is returned. 
            if delta > tol
                u1 = (-B -sqrt(delta)) / (2*A);
                u2 = (-B +sqrt(delta)) / (2*A);
                
                %want the nearest intersection, but in the direction the
                %ray is pointing.  If the nearest intersection is within
                %the tolerance, assume the ray already saw it.  
                if ((u1 > u2 || u1 < tol) && u2 > tol  )  
                    location= r.location+u2*r.direction;
                end
                if ((u1 < u2 || u2 < tol )&& u1 > tol )
                   location= r.location+u1*r.direction;
                end 
            else
                if abs(delta) < tol
                % delta around zero: find unique root.
                u = -B/(2*A);    
                location = r.location + u*r.direction;
                end
            end
            %what is the normal? 
            eCoord=(location-loc);
            if length(eCoord)==2
                normal=[eCoord(1)/b^2, eCoord(2)/a^2];
            end
            if length(eCoord)==3
                normal=[eCoord(1)/b^2, eCoord(2)/b^2, eCoord(3)/a^2];
            end
            normal=normal./sqrt(normal*normal');
             
        end
        
    end
end