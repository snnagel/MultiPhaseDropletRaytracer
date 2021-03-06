classdef TriPhaseDrop

	properties
        location
		radius
		innerRadius1
        innerRadius2
		index1
		index2
        index3
		vol1
        vol2
        vol3
		d1
        d2
        rotAngle %in degrees
	end
	
	methods
		function d=TriPhaseDrop(loc, Rd, Ri1, Ri2, i1, i2, i3, v1, v2, rotAng) %Constructor
            d.location=loc;
            d.radius=Rd;
            d.innerRadius1=Ri1;
            d.innerRadius2=Ri2;
            d.index1=i1;
            d.index2=i2;
            d.index3=i3;
            d.vol1=v1; %volume FRACTION of total drop volume must be less than 1
            d.vol2=v2;
            d.vol3=1-v1-v2
            if Ri1<0
                message='Ri1 should be > 0.  expect errors!'
            end
            %need to come back, for now assume Ri1 and Ri2 >0
            if Ri1>0 && Ri2>0 
                d.d1=fsolve(@(x)dropShapeSolver(x, Rd, Ri1, v1/(v2+d.vol3)), Ri1);
                d.d2=fsolve(@(x)dropShapeSolver(x, Rd, Ri2, (v1+v2)/d.vol3), Ri2);            
            else
                d.d1=sign(Ri1)*fsolve(@(x)dropShapeSolver(x, Rd, abs(Ri1), v1/(v2+d.vol3)), abs(Ri1));
                d.d2=sign(Ri2)*fsolve(@(x)dropShapeSolver(x, Rd, abs(Ri2), d.vol3/(v1+v2)), abs(Ri2));
            end
            if nargin>9
                d.rotAngle=rotAng;
            else
                d.rotAngle=0;
            end

        end
        
			
        function drawDrop(d, N, edgeC, c1, c2, c3)  %Draw the edges of the drop.
            if nargin<2
                N=100; %Number of points on edge of drop (increase to make smoother looking)
            end
            
            
            Rd=d.radius;
            Ri1=d.innerRadius1;
            Ri2=d.innerRadius2;
            X1=d.location(1);
            X2=d.location(2);
            %rot=d.rotAngle;
            
            
            
            intersection1=1/(2*h)*(Rd^2+d.d1^2-Ri1^2); % distance from center of drop to intersection plane of drops
            intersection2=1/(2*d.d2)*(Rd^2+d.d2^2-Ri2^2);
            
            if nargin<3
                edgeC=[0, 0, 0]
            end
            
            
            %outer color:
            if nargin<4
               c1=[1, 0.4, 0.7]
            end
            
            if nargin<5
               c2=[0.6, 0.5, 1]
            end
            
            
            if nargin<6
               c3=[1, 0.5, 1]
            end
            alphaVal=0.3
    
            thetaM_In1=acos(abs((d.d1-intersection1)/Ri));  %max angle to plot for inner phase
            thetaM_Out1=acos(intersection1/Rd) %max angle of outer phase  
            
            thetaM_In2=acos(abs((d.d2-intersection2)/Ri));  %max angle to plot for inner phase
            thetaM_Out2=acos(intersection2/Rd) %max angle of outer phase  
            
            if intersection1 > d.d1
                thetaM_In1=pi-thetaM_In1;

            end
            
            if intersection2 > d.d2
                thetaM_In2=pi-thetaM_In2;

            end            
           
            %2D Case:
            if length(d.location)==2
                
                %The outer drop: Just plot a sphere. 
                
                theta=linspace(0, 2*pi,N);
                x1=Rd*cos(theta)+X1;
                x2=Rd*sin(theta)+X2 ;
                plot(x1, x2, 'k' ,'Linewidth', 2);
                hold on
                
                %The inner drop
                %Plot a sphere up to the interface.
                if Ri1>0 && Ri2>0
                    intersection1=1/(2*d.d1)*(Rd^2+d.d1^2-Ri1^2); % distance from center of drop to intersection plane of drops
                    thetaM1=acos(abs((d.d1-intersection1)/Ri1));  %max angle to plot
                    if intersection1 > d.d1
                    thetaM1=pi-thetaM1;
                    end
                    thetai=linspace(-thetaM1, thetaM1,N);
                    x1i=-Ri1.*(cos(thetai).*sind(d.rotAngle)+sin(thetai).*cosd(d.rotAngle))+X1+d.d1.*sind(d.rotAngle);
                    x2i=Ri1.*(sin(thetai).*sind(d.rotAngle)-cos(thetai).*cosd(d.rotAngle))+X2+d.d1.*cosd(d.rotAngle);
                    plot(x1i, x2i, 'k', 'Linewidth', 2);
                    
                     % distance from center of drop to intersection plane of drops
                    thetaM2=acos(abs((d.d2-intersection2)/Ri2));  %max angle to plot
                    if intersection2 > d.d2
                    thetaM2=pi-thetaM2;
                    end
                    thetai=linspace(-thetaM2, thetaM2,N);
                    x1i=-Ri2.*(cos(thetai).*sind(d.rotAngle)+sin(thetai).*cosd(d.rotAngle))+X1+d.d2.*sind(d.rotAngle);
                    x2i=Ri2.*(sin(thetai).*sind(d.rotAngle)-cos(thetai).*cosd(d.rotAngle))+X2+d.d2.*cosd(d.rotAngle);
                    plot(x1i, x2i, 'k', 'Linewidth', 2);
                    
                              
                elseif Ri1>0 && Ri2<0
                    intersection1=1/(2*d.d1)*(Rd^2+d.d1^2-Ri1^2); % distance from center of drop to intersection plane of drops
                    thetaM1=acos(abs((d.d1-intersection1)/Ri1));  %max angle to plot
                    if intersection1 > d.d1
                    thetaM1=pi-thetaM1;
                    end
                    thetai=linspace(-thetaM1, thetaM1,N);
                    x1i=-Ri1.*(cos(thetai).*sind(d.rotAngle)+sin(thetai).*cosd(d.rotAngle))+X1+d.d1.*sind(d.rotAngle);
                    x2i=Ri1.*(sin(thetai).*sind(d.rotAngle)-cos(thetai).*cosd(d.rotAngle))+X2+d.d1.*cosd(d.rotAngle);
                    plot(x1i, x2i, 'k', 'Linewidth', 2);
                    
                    intersection2=1/(-2*d.d2)*(Rd^2+d.d2^2-Ri2^2); % distance from center of drop to intersection plane of drops
                    thetaM2=acos(abs((-d.d2-intersection2)/-Ri2));  %max angle to plot
                    if intersection2 > -d.d2
                    thetaM2=-thetaM2+pi;
                    end
                    thetai=linspace(-thetaM2, thetaM2,N);
                    x1i=Ri2.*(cos(thetai).*sind(d.rotAngle+180)+sin(thetai).*cosd(d.rotAngle+180))+X1-d.d2.*sind(d.rotAngle+180);
                    x2i=-Ri2.*(sin(thetai).*sind(d.rotAngle+180)-cos(thetai).*cosd(d.rotAngle+180))+X2-d.d2.*cosd(d.rotAngle+180);
                    plot(x1i, x2i, 'k', 'Linewidth', 2);
                    
                end
            end
            
%             %3D Case: 
%             if( length(d.location)==3);
%                 
%                 X3=d.location(3);
%                 %innerCent=d.location+[d.offCenterDist*sind(rot) ,0,d.offCenterDist*cosd(rot) ];
%                 
%                 if h>0
%                     
%                     
%                 intersection2=1/(2*h)*(Rd^2+h^2-Ri^2); % distance from center of drop to intersection plane of drops
%                     
%                 thetaMOuterPhase=acos(intersection2/Rd); %max angle of outer phase  
%                 
%                 thetaM=acos(abs((h-intersection2)/Ri));  %max angle to plot
%                     
%                 if intersection2 > h
%                     thetaM=pi-thetaM;
%                 end
%                 
%                 
%                 %The outer drop: 
%                 Phi=linspace(-pi/2,pi/2-thetaMOuterPhase, N);
%                 theta=linspace(0, 2*pi,N);
%                 [Phi,theta]=meshgrid(Phi,theta);
% 
%                 [X,Y,Z]=sph2cart(theta,Phi,Rd);
%      
%                 Xout=cosd(rot)*X-sind(rot)*Z+X1;
%                 Zout=sind(rot)*X+cosd(rot)*Z+X3;
%                 Yout=Y+X2;
%                 outerSurface=surf(Xout,Yout,Zout,'FaceColor','blue','EdgeColor','none')
%                 alpha(outerSurface, 0.15)
%                 hold on
%                 
%                 %The inner drop
%                 %Plot a sphere up to the interface. 
%                 
% 
%                     
% 
%                     thetai=linspace(-pi/2, thetaM-pi/2,N);
%                     phii=linspace(0, 2*pi, N);
%                     [phii, thetai]=meshgrid(thetai, phii);
%                     [Xtemp,Ytemp,Ztemp]=sph2cart(thetai,phii,Ri);
%                     %displace the cup:
%                     %Xtemp=Xtemp;
%                     %Ytemp=Ytemp;
%                     Ztemp=Ztemp+h;
%                     
%                     Xint=cosd(rot)*Xtemp-sind(rot)*Ztemp+X1;
%                     Zint=sind(rot)*Xtemp+cosd(rot)*Ztemp+X3;
%                     Yint=Ytemp+X2;
%                     
%                     
%                     %rotate around droplet center: 
%                     innerSurface=surf(Xint,Yint,Zint,'FaceColor',[1, 0, 0],'EdgeColor','none')
%                     alpha(innerSurface, 0.2);
%                     
%                     thetaiOuter=linspace(pi/2, thetaMOuterPhase+pi/2,N);
%                     phii=linspace(0, 2*pi, N);
%                     [phii, thetaiOuter]=meshgrid(thetaiOuter, phii);
%                     [Xtemp,Ytemp,Ztemp]=sph2cart(thetai,phii,Rd);
%                                       
%                     Xint=cosd(rot)*Xtemp-sind(rot)*Ztemp+X1;
%                     Zint=sind(rot)*Xtemp+cosd(rot)*Ztemp+X3;
%                     Yint=Ytemp+X2;
%                     
%                     innerSurface=surf(Xint,Yint,Zint,'FaceColor','red','EdgeColor','none')
%                     alpha(innerSurface, 0.2);
%                     
% 
%                 else
%                     message= 'CODE UNFINISHED DON"T TRUST PLOT'
%                     intersection2=-1/(2*h)*(Rd^2+h^2-Ri^2); % distance from center of drop to intersection plane of drops
%                     
%                     thetaM=acos((-h-intersection2)/Ri);
%                     if intersection2 > h
%                     thetaM=thetaM-pi/2;
%                     end
%                     thetai=linspace(-pi/2, thetaM-pi/2,N);
%                     phii=linspace(0, 2*pi, N);
%                     [phii, thetai]=meshgrid(thetai, phii);
%                     [Xtemp,Ytemp,Ztemp]=sph2cart(thetai,phii,Ri);
%                     %displace the cup:
%                     Xtemp=Xtemp;
%                     Ytemp=Ytemp;
%                     Ztemp=Ztemp+h;
%                     
%                     Xint=cosd(rot)*Xtemp-sind(rot)*Ztemp+X1;
%                     Zint=sind(rot)*Xtemp+cosd(rot)*Ztemp+X3;
%                     Yint=Ytemp+X2;
%                     
%                     
%                     %rotate around droplet center: 
%                     innerSurface=surf(Xint,Yint,Zint,'FaceColor','blue','EdgeColor','none');
%                     alpha(innerSurface, 0.4);
%                 end
%     
% 
% 
% 
%
%            end
        end
        
        
        
        function drawDropB(d, N)  %Draw the edges of the drop.
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
                
                %The outer drop: Just plot a sphere. 
                
                theta=linspace(0, 2*pi,N);
                x1=Rd*cos(theta)+X1;
                x2=Rd*sin(theta)+X2 ;
                plot(x1, x2, 'w' ,'Linewidth', 2);
                hold on
                
                %The inner drop
                %Plot a sphere up to the interface. 
                if h>0
                    intersection=1/(2*h)*(Rd^2+h^2-Ri^2); % distance from center of drop to intersection plane of drops
                    thetaM=acos(abs((h-intersection)/Ri));  %max angle to plot
                    if intersection > h
                    thetaM=pi-thetaM;
                    end
                    thetai=linspace(-thetaM, thetaM,N);
                    x1i=-Ri.*(cos(thetai).*sind(d.rotAngle)+sin(thetai).*cosd(d.rotAngle))+X1+h.*sind(d.rotAngle);
                    x2i=Ri.*(sin(thetai).*sind(d.rotAngle)-cos(thetai).*cosd(d.rotAngle))+X2+h.*cosd(d.rotAngle);
                    plot(x1i, x2i, 'w', 'Linewidth', 2);
                else
                    message='Choose Positive Internal Curvature'
%                     intersection=-1/(2*h)*(Rd^2+h^2-Ri^2); % distance from center of drop to intersection plane of drops
%                     
%                     thetaM=acos((-h-intersection)/Ri)
%                     if h > intersection 
%                     thetaM=pi-thetaM;
%                     end
%                     thetai=linspace(-thetaM, thetaM,N);
%                     x1i=Ri*(cos(thetai).*cosd(d.rotAngle)+sin(thetai).*sind(d.rotAngle))+X1+h.*cosd(d.rotAngle);
%                     x2i=Ri*(sin(thetai).*cosd(d.rotAngle)-cos(thetai).*sind(d.rotAngle))+X2+h.*sind(d.rotAngle);
%                     plot(x1i, x2i, 'Linewidth', 2);
                end
                else
                
                    
            
            end
        end
        
        
        
%         function drawDrop(d, N)  %Draw the edges of the drop.
%             if nargin<2
%                 N=100; %Number of points on edge of drop (increase to make smoother looking)
%             end
%             Rd=d.radius;
%             Ri=d.innerRadius;
%             h=d.offCenterDist
%             X1=d.location(1);
%             X2=d.location(2);
%             rot=d.rotAngle;
%             %2D Case:
%             if length(d.location)==2
%                 
%                 %The outer drop: Just plot a sphere. 
%                 
%                 theta=linspace(0, 2*pi,N);
%                 x1=Rd*cos(theta)+X1;
%                 x2=Rd*sin(theta)+X2 ;
%                 plot(x1, x2, 'k' ,'Linewidth', 2);
%                 hold on
%                 
%                 %The inner drop
%                 %Plot a sphere up to the interface. 
%                 if h>0
%                     intersection=1/(2*h)*(Rd^2+h^2-Ri^2); % distance from center of drop to intersection plane of drops
%                     thetaM=acos(abs((h-intersection)/Ri));  %max angle to plot
%                     if intersection > h
%                     thetaM=pi-thetaM;
%                     end
%                     thetai=linspace(-thetaM, thetaM,N);
%                     x1i=-Ri.*(cos(thetai).*sind(d.rotAngle)+sin(thetai).*cosd(d.rotAngle))+X1+h.*sind(d.rotAngle);
%                     x2i=Ri.*(sin(thetai).*sind(d.rotAngle)-cos(thetai).*cosd(d.rotAngle))+X2+h.*cosd(d.rotAngle);
%                     plot(x1i, x2i, 'k', 'Linewidth', 2);
%                 else
%                     message='Choose Positive Internal Curvature'
% %                     intersection=-1/(2*h)*(Rd^2+h^2-Ri^2); % distance from center of drop to intersection plane of drops
% %                     
% %                     thetaM=acos((-h-intersection)/Ri)
% %                     if h > intersection 
% %                     thetaM=pi-thetaM;
% %                     end
% %                     thetai=linspace(-thetaM, thetaM,N);
% %                     x1i=Ri*(cos(thetai).*cosd(d.rotAngle)+sin(thetai).*sind(d.rotAngle))+X1+h.*cosd(d.rotAngle);
% %                     x2i=Ri*(sin(thetai).*cosd(d.rotAngle)-cos(thetai).*sind(d.rotAngle))+X2+h.*sind(d.rotAngle);
% %                     plot(x1i, x2i, 'Linewidth', 2);
%                 end
%                 else
%                 
%                     
%             
%             end
%             
            %3D Case: 
%             if( length(d.location)==3);
%                 
%                 X3=d.location(3);
%                 %innerCent=d.location+[d.offCenterDist*sind(rot) ,0,d.offCenterDist*cosd(rot) ];
%                 
%                 if h>0
%                     
%                     
%                 intersection=1/(2*h)*(Rd^2+h^2-Ri^2); % distance from center of drop to intersection plane of drops
%                     
%                 thetaMOuterPhase=acos(intersection/Rd); %max angle of outer phase  
%                 
%                 thetaM=acos(abs((h-intersection)/Ri));  %max angle to plot
%                     
%                 if intersection > h
%                     thetaM=pi-thetaM;
%                 end
%                 
%                 
%                 %The outer drop: 
%                 Phi=linspace(-pi/2,pi/2-thetaMOuterPhase, N);
%                 theta=linspace(0, 2*pi,N);
%                 [Phi,theta]=meshgrid(Phi,theta);
% 
%                 [X,Y,Z]=sph2cart(theta,Phi,Rd);
%      
%                 Xout=cosd(rot)*X-sind(rot)*Z+X1;
%                 Zout=sind(rot)*X+cosd(rot)*Z+X3;
%                 Yout=Y+X2;
%                 outerSurface=surf(Xout,Yout,Zout,'FaceColor','blue','EdgeColor','none')
%                 alpha(outerSurface, 0.15)
%                 hold on
%                 
%                 %The inner drop
%                 %Plot a sphere up to the interface. 
%                 
% 
%                     
% 
%                     thetai=linspace(-pi/2, thetaM-pi/2,N);
%                     phii=linspace(0, 2*pi, N);
%                     [phii, thetai]=meshgrid(thetai, phii);
%                     [Xtemp,Ytemp,Ztemp]=sph2cart(thetai,phii,Ri);
%                     %displace the cup:
%                     %Xtemp=Xtemp;
%                     %Ytemp=Ytemp;
%                     Ztemp=Ztemp+h;
%                     
%                     Xint=cosd(rot)*Xtemp-sind(rot)*Ztemp+X1;
%                     Zint=sind(rot)*Xtemp+cosd(rot)*Ztemp+X3;
%                     Yint=Ytemp+X2;
%                     
%                     
%                     %rotate around droplet center: 
%                     innerSurface=surf(Xint,Yint,Zint,'FaceColor',[1, 0, 0],'EdgeColor','none')
%                     alpha(innerSurface, 0.2);
%                     
%                     thetaiOuter=linspace(pi/2, thetaMOuterPhase+pi/2,N);
%                     phii=linspace(0, 2*pi, N);
%                     [phii, thetaiOuter]=meshgrid(thetaiOuter, phii);
%                     [Xtemp,Ytemp,Ztemp]=sph2cart(thetai,phii,Rd);
%                                       
%                     Xint=cosd(rot)*Xtemp-sind(rot)*Ztemp+X1;
%                     Zint=sind(rot)*Xtemp+cosd(rot)*Ztemp+X3;
%                     Yint=Ytemp+X2;
%                     
%                     innerSurface=surf(Xint,Yint,Zint,'FaceColor','red','EdgeColor','none')
%                     alpha(innerSurface, 0.2);
%                     
% 
%                 else
%                     message= 'CODE UNFINISHED DON"T TRUST PLOT'
%                     intersection=-1/(2*h)*(Rd^2+h^2-Ri^2); % distance from center of drop to intersection plane of drops
%                     
%                     thetaM=acos((-h-intersection)/Ri);
%                     if intersection > h
%                     thetaM=thetaM-pi/2;
%                     end
%                     thetai=linspace(-pi/2, thetaM-pi/2,N);
%                     phii=linspace(0, 2*pi, N);
%                     [phii, thetai]=meshgrid(thetai, phii);
%                     [Xtemp,Ytemp,Ztemp]=sph2cart(thetai,phii,Ri);
%                     %displace the cup:
%                     Xtemp=Xtemp;
%                     Ytemp=Ytemp;
%                     Ztemp=Ztemp+h;
%                     
%                     Xint=cosd(rot)*Xtemp-sind(rot)*Ztemp+X1;
%                     Zint=sind(rot)*Xtemp+cosd(rot)*Ztemp+X3;
%                     Yint=Ytemp+X2;
%                     
%                     
%                     %rotate around droplet center: 
%                     innerSurface=surf(Xint,Yint,Zint,'FaceColor','blue','EdgeColor','none');
%                     alpha(innerSurface, 0.4);
%                 end
%     
% 
% 
% 
% 
%             end
%         end
%         
        
        
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
                if h>0
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
                else
                
                    
            
            end
        end
        
        
        
        %returns the tangent plane from the drop, in the direction dir. The
        %plane is in the form %p(2)*x1+p(3)*x2+p(4)*x3= p(1)
        function p=tangentPlane(dir, d)
            dir=dir./sqrt(dir*dir'); %make sure its normalized. 
            p=[d.radius+dir*d.location', dir];
            
        end 


    end
end

		