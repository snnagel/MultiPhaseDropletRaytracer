classdef dropElliptical

	properties
        location
		radius
		a
        b
		innerIndex
		outerIndex
		volumeRatio
		offCenterDist
        rotAngle %in degrees
        interfaceZ  %z coordinate of the three phase contact line (relative to center of drop). 
	end
	
	methods
		function d=dropElliptical(loc, Rd, a, b, outerI, innerI, vr, rotAng) %Constructor
            d.location=loc;
            d.radius=Rd;
            d.outerIndex=outerI;
            d.innerIndex=innerI;
            d.volumeRatio=vr;

            d.a=a;
            d.b=b;
            
            %Determine start values using spherical case: 
            if a > b
                h0=fsolve(@(x)dropShapeSolver(x, Rd, a, vr), a);
            else   
            h0=fsolve(@(x)dropShapeSolver(x, Rd, b, vr), b);
            end
            l0=1/(2*h0)*(Rd^2+h0^2-a^2);
            %Use those values as start point to determine elliptical case:
           
            %solveOptions = optimoptions(@fsolve, 'MaxFunctionEvaluations', 1000, 'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'Algorithm', 'Levenberg-Marquardt')
            shape=fsolve(@(x)dropShapeSolverElliptical(x, Rd, a, b, vr), [h0, l0]);

            d.offCenterDist=shape(1);
            d.interfaceZ=shape(2);
            
            
            if nargin>7
                d.rotAngle=rotAng;
            else
                d.rotAngle=0;
            end

                
        end
        
			
		function n= index(d, loc) %returns index at that location assuming one drop present
            'Sorry, function index has not been implemented yet' 
            n=nan;
        end
        function drawDropE(d)  %Draw the edges of the drop.
            N=400; %Number of points on edge of drop (increase to make smoother looking)
            Rd=d.radius;
            a=d.a;
            b=d.b;
            h=d.offCenterDist;
            X1=d.location(1);
            X2=d.location(2);
            rot=d.rotAngle;
            %2D Case:
            if length(d.location)==2
                N=4000
                %The outer drop: Just plot a sphere. 
                
                theta=linspace(0, 2*pi,N);
                x1=Rd*cos(theta)+X1;
                x2=Rd*sin(theta)+X2 ;
                plot(x1, x2, 'Linewidth', 2);
                hold on
                
                %The inner drop
                %Plot a sphere up to the interface. 
                if h>0
                    intersection=((b^2*h - sqrt((-Rd^2 + h^2)*a^2*b^2 + a^2*b^4 + Rd^2*a^4 - a^4*b^2))/(b^2 - a^2)); % distance from center of drop to intersection plane of drops
                    thetaM=acos(abs((h-intersection)/a));  %max angle to plot
                    if intersection > h
                    thetaM=pi-thetaM;
                    end
                    thetai=linspace(-thetaM, thetaM,N);
                    x1i=b*sin(thetai);
                    x2i=-a*cos(thetai)+h;
                    plot(x1i, x2i, 'Linewidth', 2);
                end
                
                    
            
            end
            
            %3D Case: 
            if( length(d.location)==3);
                N=300
                X3=d.location(3);
                %innerCent=d.location+[d.offCenterDist*sind(rot) ,0,d.offCenterDist*cosd(rot) ];
                
                if h>0
                    
                    
                intersection=((b^2*h - sqrt((-Rd^2 + h^2)*a^2*b^2 + a^2*b^4 + Rd^2*a^4 - a^4*b^2))/(b^2 - a^2)); % distance from center of drop to intersection plane of drops
                    
                thetaMOuterPhase=acos(intersection/Rd); %max angle of outer phase  
                
                thetaM=acos(abs((h-intersection)/b));  %max angle to plot
                    
                if intersection > h
                    thetaM=thetaM-pi/2;
                end
                
                
                %The outer drop: 
                Phi=linspace(-pi/2,pi/2-thetaMOuterPhase, N);
                theta=linspace(0, 2*pi,N);
                [Phi,theta]=meshgrid(Phi,theta);

                [X,Y,Z]=sph2cart(theta,Phi,Rd);
     
                Xout=cosd(rot)*X-sind(rot)*Z+X1;
                Zout=sind(rot)*X+cosd(rot)*Z+X3;
                Yout=Y+X2;
                outerSurface=surf(Xout,Yout,Zout,'FaceColor','blue','EdgeColor','none')
                alpha(outerSurface, 0.15)
                hold on
                
                %The inner drop
                %Plot a sphere up to the interface. 
                

                    

                    thetai=linspace(-pi/2, thetaM-pi/2,N);
                    phii=linspace(0, 2*pi, N);
                    [thetai, phii]=meshgrid(thetai, phii);
                    ri=sqrt(1./(sin(thetai).^2./a^2+cos(thetai).^2./b^2));
                    %ri=a;
                    [Xtemp,Ytemp,Ztemp]=sph2cart(phii,thetai,ri);
                    %displace the cup:
                    %Xtemp=Xtemp;
                    %Ytemp=Ytemp;
                    %Ztemp=Ztemp+h;
                    
                    Xint=Xtemp+X1;
                    Zint=Ztemp+h+X3;
                    Yint=Ytemp+X2;
                    
                    
                    %rotate around droplet center: 
                    innerSurface=surf(Xint,Yint,Zint,'FaceColor',[1, 0, 0],'EdgeColor','none')
                    alpha(innerSurface, 0.2)
                    
                    thetaiOuter=linspace(pi/2, thetaMOuterPhase+pi/2,N);
                    phii=linspace(0, 2*pi, N);
                    [phii, thetaiOuter]=meshgrid(thetaiOuter, phii);
                    [Xtemp,Ytemp,Ztemp]=sph2cart(thetaiOuter,phii,Rd);
                                      
                    Xint=cosd(rot)*Xtemp-sind(rot)*Ztemp+X1;
                    Zint=sind(rot)*Xtemp+cosd(rot)*Ztemp+X3;
                    Yint=Ytemp+X2;
                    
                    innerSurface=surf(Xint,Yint,Zint,'FaceColor','red','EdgeColor','none')
                    alpha(innerSurface, 0.2)
                    

                else
                    message= 'CODE UNFINISHED DON"T TRUST PLOT'
                    intersection=-1/(2*h)*(Rd^2+h^2-Ri^2); % distance from center of drop to intersection plane of drops
                    
                    thetaM=acos((-h-intersection)/Ri)
                    if intersection > h
                    thetaM=thetaM-pi/2;
                    end
                    thetai=linspace(-pi/2, thetaM-pi/2,N);
                    phii=linspace(0, 2*pi, N);
                    [phii, thetai]=meshgrid(thetai, phii);
                    [Xtemp,Ytemp,Ztemp]=sph2cart(thetai,phii,Ri);
                    %displace the cup:
                    Xtemp=Xtemp;
                    Ytemp=Ytemp;
                    Ztemp=Ztemp+h;
                    
                    Xint=cosd(rot)*Xtemp-sind(rot)*Ztemp+X1;
                    Zint=sind(rot)*Xtemp+cosd(rot)*Ztemp+X3;
                    Yint=Ytemp+X2;
                    
                    
                    %rotate around droplet center: 
                    innerSurface=surf(Xint,Yint,Zint,'FaceColor','blue','EdgeColor','none')
                    alpha(innerSurface, 0.4)
                end
               
            

            
            
 

            camlight left; 
            lighting phong



            end
        end
        

        function relativeCoord=ThreePhaseContactLineEl(d)
            %Returns the location relative to the center of the drop that the
            %interface occurs in local droplet coordinates.
            %Note in 2D and 3D value returned is x, z location of interface
            %where z is axis of droplet, and x is distance from that axis. The
            %origin is at the center of the droplet. 
            Rd=d.radius;
            %WARNING UNTESTED
            h=d.offCenterDist;
            Zi=(d.b^2*h-sqrt((-Rd^2 + h^2)*d.a^2*d.b^2 + d.a^2*d.b^4 + Rd^2*d.a^4 - d.a^4*d.b^2))/(d.b^2 - d.a^2);
            Xi=sqrt(Rd^2-Zi^2);
            relativeCoord=[ Zi, Xi];        
        end
            
    end
end

		