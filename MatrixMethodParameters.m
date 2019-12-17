%Focal lenght as function of curvature (NOT contact angle)
%paraxial approximation
%ray transfer method
clear all
close all
%Indexes of refraction: 
nm=1.33; %surrounding Medium-water
n1=1.39;%innerPhase
n2=1.27; %outerPhase
angle=0
figure

%DropRadii=[1.1:.1:4, 4.5:.5:20]
Rd=1
P=[]
fParaxVr=[]
volRatios=[.5:.1:2]
%volRatios=[0.5:.1:1, 1.2:.2:2]
numPoints=length(volRatios)

contactAngle=pi/4

P=[]

fParaxVr=nan(numPoints, 1)
for n1=[1.4, 1.5, 1.6]
    subplot(2, 2, 2)
for m=[1:numPoints]
    vr=volRatios(m)
    

[solu,fval,exitflag,output]=fsolve(@(x)dropShapeSolverContactAngle(x, Rd, vr, cos(contactAngle)), [0.8,1] )
if (exitflag<0)
    continue
end
d=solu(1)
Ri=solu(2)
theDrop=drop([Rd+1,0], Rd, Ri, n1, n2, vr, angle);

% figure
% drawDrop(theDrop)
% axis image
% title(['contactAngle=', num2str(contactAngle), '  and vr =', num2str(vr)])

l1 = Rd - d + Ri;
l2 = Rd - Ri + d;

intoDrop = [[1,  0]; [(nm - n1)/(Rd*n1), nm/n1]];
throughFirstPhase = [[1, l1]; [0, 1]];
intoInterface = [[1,  0]; [(n1 - n2)/(-Ri*n2), n1/n2]];
throughSecondPhase = [[1, l2]; [0, 1]];
outOfDrop = [[1,  0]; [(n2 - nm)/(-1*Rd*nm), n2/nm]];

dropMat=outOfDrop*throughSecondPhase*intoInterface*throughFirstPhase*intoDrop


s=-dropMat(1,1)/dropMat(2,1)
effFocalLength=-1/dropMat(2,1)

fParaxVr(m)=effFocalLength  %effective Focal Length
end

plot(volRatios, fParaxVr, 'Linewidth', 2)
hold on
set(0,'defaulttextinterpreter','latex')
xlabel('Volume Ratio')
ylabel('Effective Focal Length')
%print(gcf, 'varyingVolumeRatio/FocalLengthVsVolRatio', '-deps')
xlim([0.5,2])
%%
subplot(2, 2, 1)
vr=1
CAs=linspace(0.1, pi/2, numPoints)
for m=[1:numPoints]
    %vr=volRatios(m)
    contactAngle=CAs(m)

[solu,fval,exitflag,output]=fsolve(@(x)dropShapeSolverContactAngle(x, Rd, vr, cos(contactAngle)), [0.8,1] )
if (exitflag<0)
    continue
end
d=solu(1)
Ri=solu(2)
theDrop=drop([Rd+1,0], Rd, Ri, n1, n2, vr, angle);

% figure
% drawDrop(theDrop)
% axis image
% title(['contactAngle=', num2str(contactAngle), '  and vr =', num2str(vr)])

l1 = Rd - d + Ri;
l2 = Rd - Ri + d;

intoDrop = [[1,  0]; [(nm - n1)/(Rd*n1), nm/n1]];
throughFirstPhase = [[1, l1]; [0, 1]];
intoInterface = [[1,  0]; [(n1 - n2)/(-Ri*n2), n1/n2]];
throughSecondPhase = [[1, l2]; [0, 1]];
outOfDrop = [[1,  0]; [(n2 - nm)/(-1*Rd*nm), n2/nm]];

dropMat=outOfDrop*throughSecondPhase*intoInterface*throughFirstPhase*intoDrop


s=-dropMat(1,1)/dropMat(2,1)
effFocalLength=-1/dropMat(2,1)

fParaxCA(m)=effFocalLength  %effective Focal Length
end

plot(CAs*180/pi, fParaxCA, 'Linewidth', 2)
hold on
%fontsize=18
set(0,'defaulttextinterpreter','latex')
xlabel('Contact Angle')
ylabel('Effective Focal Length')
xlim([8, 90])
%legend( '\theta_c =0.1', '\theta_c =\pi/6', '\theta_c =\pi/4')
%print(gcf, 'varyingVolumeRatio/FocalLengthVsVolRatio', '-deps')
%xlim([0.5,2])

%%
subplot(2, 2, 4)
vr=1
contactAngle=pi/4
n2s=linspace(1.2, 1.6, numPoints)
for m=[1:numPoints]
    %vr=volRatios(m)
    %contactAngle=CAs(m)
    n2=n2s(m)

[solu,fval,exitflag,output]=fsolve(@(x)dropShapeSolverContactAngle(x, Rd, vr, cos(contactAngle)), [0.8,1] )
if (exitflag<0)
    continue
end
d=solu(1)
Ri=solu(2)
theDrop=drop([Rd+1,0], Rd, Ri, n1, n2, vr, angle);

% figure
% drawDrop(theDrop)
% axis image
% title(['contactAngle=', num2str(contactAngle), '  and vr =', num2str(vr)])

l1 = Rd - d + Ri;
l2 = Rd - Ri + d;

intoDrop = [[1,  0]; [(nm - n1)/(Rd*n1), nm/n1]];
throughFirstPhase = [[1, l1]; [0, 1]];
intoInterface = [[1,  0]; [(n1 - n2)/(-Ri*n2), n1/n2]];
throughSecondPhase = [[1, l2]; [0, 1]];
outOfDrop = [[1,  0]; [(n2 - nm)/(-1*Rd*nm), n2/nm]];

dropMat=outOfDrop*throughSecondPhase*intoInterface*throughFirstPhase*intoDrop


s=-dropMat(1,1)/dropMat(2,1)
effFocalLength=-1/dropMat(2,1)

fParaxN2(m)=effFocalLength  %effective Focal Length
end
%%
plot(n2s, fParaxN2, 'Linewidth', 2)
hold on
%fontsize=18
xlabel('$n_2$')
ylabel('Effective Focal Length')
%legend( '\theta_c =0.1', '\theta_c =\pi/6', '\theta_c =\pi/4')
%print(gcf, 'varyingVolumeRatio/FocalLengthVsVolRatio', '-deps')
%xlim([0.5,2])
ylim([-5, 20])


end
subplot(2, 2, 1)
legend('n_1=1.4','n_1=1.5','n_1=1.6', 'location', 'best')
ylim([-5, 20])
subplot(2, 2, 2)
legend('n_1=1.4','n_1=1.5','n_1=1.6', 'location', 'best')
ylim([-5, 20])
subplot(2, 2, 4)
legend('n_1=1.4','n_1=1.5','n_1=1.6', 'location', 'best')
ylim([-5, 20])
%%
subplot(2, 2, 3)
kk=1
for n2=[1.27, 1.5]
vr=1
contactAngle=pi/4
n1s=linspace(1.33, 1.6, numPoints*100)
[solu,fval,exitflag,output]=fsolve(@(x)dropShapeSolverContactAngle(x, Rd, vr, cos(contactAngle)), [0.8,1] );
d=solu(1);
Ri=solu(2);
for m=[1:numPoints*100]
    %vr=volRatios(m)
    %contactAngle=CAs(m)
    n1=n1s(m);


theDrop=drop([Rd+1,0], Rd, Ri, n1, n2, vr, angle);

% figure
% drawDrop(theDrop)
% axis image
% title(['contactAngle=', num2str(contactAngle), '  and vr =', num2str(vr)])

l1 = Rd - d + Ri;
l2 = Rd - Ri + d;

intoDrop = [[1,  0]; [(nm - n1)/(Rd*n1), nm/n1]];
throughFirstPhase = [[1, l1]; [0, 1]];
intoInterface = [[1,  0]; [(n1 - n2)/(-Ri*n2), n1/n2]];
throughSecondPhase = [[1, l2]; [0, 1]];
outOfDrop = [[1,  0]; [(n2 - nm)/(-1*Rd*nm), n2/nm]];

dropMat=outOfDrop*throughSecondPhase*intoInterface*throughFirstPhase*intoDrop;


s=-dropMat(1,1)/dropMat(2,1);
effFocalLength=-1/dropMat(2,1);

fParaxN1(m)=effFocalLength;  %effective Focal Length

end
%%
if kk==1
plot(n1s, fParaxN1, 'Linewidth', 2, 'color', [0, 0.5, 0.5])
elseif kk==2
    plot(n1s, fParaxN1, 'Linewidth', 2, 'color', [0.5, 0, 0.5])
end
kk=kk+1
hold on
%fontsize=18
xlabel('$n_1$')
ylabel('Effective Focal Length')
%legend( '\theta_c =0.1', '\theta_c =\pi/6', '\theta_c =\pi/4')
%print(gcf, 'varyingVolumeRatio/FocalLengthVsVolRatio', '-deps')
%xlim([0.5,2])
%ylim([-5, 20])

%%
end
legend('n_2=1.27', 'n_2=1.5')
ylim([0, 200])
xlim([1.32, 1.4])
set(gcf, 'color', 'w')
