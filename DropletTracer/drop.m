classdef drop

	properties
        location %droplet center
		radius %Rd
		innerRadius %Ri
		innerIndex
		outerIndex
		volumeRatio
		offCenterDist
        rotAngle %in degrees.  Only one angle is implemented even in 3D. 
	end
	
	methods
		function d=drop(loc, Rd, Ri, outerI, innerI, vr, rotAng) %Constructor
            d.location=loc;
            d.radius=Rd;
            d.outerIndex=outerI;
            d.innerIndex=innerI;
            d.volumeRatio=vr;
            

            if Ri<0 %flip the droplet to accomodate this.  
                d.innerRadius=-Ri;
                d.offCenterDist=fsolve(@(x)dropShapeSolver(x, Rd, -1*Ri, vr), -1*Ri);
                d.outerIndex=innerI;
                d.innerIndex=outerI;
                if nargin>6
                    rotAng=rotAng+180
                else
                    rotAng=180
                end
            else
                d.innerRadius=Ri;
                d.offCenterDist=fsolve(@(x)dropShapeSolver(x, Rd, Ri, vr), Ri);
                
                if nargin>6
                    rotAng=rotAng
                else
                    rotAng=0
                end
                
            end

            d.rotAngle=rotAng;
   
        end
        		
        function drawDrop(d, N, edgeC, colorO, colorI)  %Draw the edges of the drop.
            if nargin<2
                N=100; %Number of points on edge of drop (increase to make smoother looking)
            end
            
            %default colors- higher index pink, lower purple
            
            %edge color:
            if nargin<3
                edgeC=[0, 0, 0]
            end
            
            
            %outer color:
            if nargin<4
                if d.innerIndex>d.outerIndex
                    colorO=[0.6, 0.5, 1] %fluoro
                else
                     colorO=[1, 0.4, 0.7] %hyrdo
                end
            end
            
            %inner color:
            if nargin<5
                if d.innerIndex<d.outerIndex
                    colorI=[0.6, 0.5, 1] %fluoro
                else
                    colorI=[1, 0.4, 0.7] %hyrdo
                end
            end
            

            
            alphaVal=0.3 %how transparent in 3D case.  
            
            Rd=d.radius;
            Ri=d.innerRadius;
            h=d.offCenterDist
            X1=d.location(1);
            X2=d.location(2);
            rot=d.rotAngle;
                        
            
            intersection=1/(2*h)*(Rd^2+h^2-Ri^2); % distance from center of drop to intersection plane of drops
            thetaM_In=acos(abs((h-intersection)/Ri));  %max angle to plot for inner phase
            thetaM_Out=acos(intersection/Rd) %max angle of outer phase  
            
            if intersection > h
                thetaM_In=pi-thetaM_In;
                %thetaM_Out=pi-thetaM_Out
            end
      
            %2D Case:
            if length(d.location)==2
                
                %The outer phase: 
                thetao=linspace(thetaM_Out,2*pi-thetaM_Out, N);                             
                x1O=-Rd.*(cos(thetao).*sind(d.rotAngle+180)+sin(thetao).*cosd(d.rotAngle+180))+X1;
                x2O=Rd.*(sin(thetao).*sind(d.rotAngle+180)-cos(thetao).*cosd(d.rotAngle+180))+X2;
               
                %internal interface:
                thetai=linspace(-thetaM_In, thetaM_In,N);
                x1io=-Ri.*(cos(thetai).*sind(d.rotAngle)+sin(thetai).*cosd(d.rotAngle))+X1+h.*sind(d.rotAngle);
                x2io=Ri.*(sin(thetai).*sind(d.rotAngle)-cos(thetai).*cosd(d.rotAngle))+X2+h.*cosd(d.rotAngle);
                
                %The inner phase
                thetaiO=linspace(-thetaM_Out, thetaM_Out,N);
                x1I=-Rd.*(cos(thetaiO).*sind(d.rotAngle+180)+sin(thetaiO).*cosd(d.rotAngle+180))+X1;
                x2I=Rd.*(sin(thetaiO).*sind(d.rotAngle+180)-cos(thetaiO).*cosd(d.rotAngle+180))+X2;
                        
                %plot: 
                fill([x1io, fliplr(x1O) ], [x2io, fliplr(x2O) ], colorO,  'linewidth', 2, 'edgeColor', edgeC)
                hold on
                f=fill([x1io, x1I],[x2io, x2I], colorI, 'linewidth', 2, 'edgeColor', edgeC)
                
            end
    
            
            
            
            %3D Case: 
            if( length(d.location)==3);
                
                X3=d.location(3);
                                
                %The outer drop: 
                Phi=linspace(-pi/2,pi/2-thetaM_Out, N);
                theta=linspace(0, 2*pi,N);
                [Phi,theta]=meshgrid(Phi,theta);

                [X,Y,Z]=sph2cart(theta,Phi,Rd);
     
                Xout=cosd(rot)*X-sind(rot)*Z+X1;
                Zout=sind(rot)*X+cosd(rot)*Z+X3;
                Yout=Y+X2;
                outerSurface=surf(Xout,Yout,Zout,'FaceColor',colorO,'EdgeColor','none')
                alpha(outerSurface, alphaVal)
                hold on
                
                %The inner drop

                thetai=linspace(-pi/2, thetaM_In-pi/2,N);
                phii=linspace(0, 2*pi, N);
                [phii, thetai]=meshgrid(thetai, phii);
                [Xtemp,Ytemp,Ztemp]=sph2cart(thetai,phii,Ri);
                
                %displace the cup:
                Ztemp=Ztemp+h;
                
                %rotate around droplet center: 
                Xint=cosd(rot)*Xtemp-sind(rot)*Ztemp+X1;
                Zint=sind(rot)*Xtemp+cosd(rot)*Ztemp+X3;
                Yint=Ytemp+X2;

                innerSurface=surf(Xint,Yint,Zint,'FaceColor',colorI,'EdgeColor','none')
                alpha(innerSurface, alphaVal);

                thetaiOuter=linspace(pi/2, thetaM_Out+pi/2,N);
                phii=linspace(0, 2*pi, N);
                [phii, thetaiOuter]=meshgrid(thetaiOuter, phii);
                [Xtemp,Ytemp,Ztemp]=sph2cart(thetai,phii,Rd);

                Xint=cosd(rot)*Xtemp-sind(rot)*Ztemp+X1;
                Zint=sind(rot)*Xtemp+cosd(rot)*Ztemp+X3;
                Yint=Ytemp+X2;

                innerSurface=surf(Xint,Yint,Zint,'FaceColor',colorI,'EdgeColor','none')
                alpha(innerSurface, alphaVal);

            end
        end
        
        
        
       %2D only:
        function drawDropIndexMatch(d, N)  %Draw the edges of the drop but dashed on upper interface.
            if nargin<2
                N=100; %Number of points on edge of drop (increase to make smoother looking)
            end
            Rd=d.radius;
            Ri=d.innerRadius;
            h=d.offCenterDist
            X1=d.location(1);
            X2=d.location(2);
            rot=d.rotAngle;
            %2D Case:
            if length(d.location)==2  
                %The inner drop
                %Plot a sphere up to the interface. 
                    intersection=1/(2*h)*(Rd^2+h^2-Ri^2); % distance from center of drop to intersection plane of drops
                    thetaM=acos(abs((h-intersection)/Ri));  %max angle to plot
                    if intersection > h
                    thetaM=pi-thetaM;
                    end
                    thetai=linspace(-thetaM, thetaM,N);
                    x1i=-Ri.*(cos(thetai).*sind(d.rotAngle)+sin(thetai).*cosd(d.rotAngle))+X1+h.*sind(d.rotAngle);
                    x2i=Ri.*(sin(thetai).*sind(d.rotAngle)-cos(thetai).*cosd(d.rotAngle))+X2+h.*cosd(d.rotAngle);
                    plot(x1i, x2i, 'w', 'Linewidth', 2);
                    hold on
                    %outer interfaces:
                    intersectionO=1/(2*h)*(Ri^2+h^2-Rd^2);
                    thetaIND=acos(abs((h-intersectionO)/Rd));  %max angle to plot
                    theta=linspace(-thetaIND, thetaIND,N);
                    x1=Rd*sin(theta)+X1;
                    x2=Rd*cos(theta)+X2 ;
                    plot(x1, x2, 'w:' ,'Linewidth', 2);
                    
                    %outer:
                    theta=linspace(thetaIND, 2*pi-thetaIND,N);
                    x1=Rd*sin(theta)+X1;
                    x2=Rd*cos(theta)+X2 ;
                    plot(x1, x2, 'w' ,'Linewidth', 2);
            end
        end

        %returns the tangent plane from the drop, in the direction dir. The
        %plane is in the form %p(2)*x1+p(3)*x2+p(4)*x3= p(1)
        function p=tangentPlane(dir, d)
            dir=dir./sqrt(dir*dir'); %make sure its normalized. 
            p=[d.radius+dir*d.location', dir];
            
        end 

        function relativeCoord=ThreePhaseContactLine(d)
            %Returns the location relative to the center of the drop that the
            %interface occurs in local droplet coordinates.
            %Note in 2D and 3D value returned is x, z location of interface
            %where z is axis of droplet, and x is distance from that axis. The
            %origin is at the center of the droplet. 
            Rd=d.radius;
            Ri=d.innerRadius;
            h=d.offCenterDist;
            Xi=1/(2*h)*(Rd^2+h^2-Ri^2);
            Zi=sqrt(Rd^2-Xi^2);
            relativeCoord=[Xi, Zi];        
        end
    end
end

		