close all
addpath('C:\Users\Sara\Dropbox (MIT)\DropletRayTracer')

%Hexagon: 
N=5
n1=1.51
n2=1
figure
D=20
DrawPoly(N, D);
x0=-0.85*D/2;
thetaIn=5 %degrees

for x0=[-D/2:.5:0]
testRay=ray([x0, 1], [sind(thetaIn), -cosd(thetaIn)], n1, 1, [0, 0, 1]);

testOut=tracePoly(testRay, N, D, n2);
if (testOut.amplitude>0.1)
drawRay(propogate(testOut, 5));
end

end
%%
N=3
D=20
n1=1.151
n2=1
close all
figure
for ii=1:3
for jj=0:5
    N=ii+2
    subplot(4, 6, (ii-1)*6+jj+1)
    DrawPoly(N, D);
    title(['\theta_i=' num2str(jj) '^o'])
x0=-0.85*D/2;
thetaIn=jj %degrees

for x0=[-D/2:.5:0]
testRay=ray([x0, 1], [sind(thetaIn), -cosd(thetaIn)], n1, 1, [0, 0, 1]);

testOut=tracePoly(testRay, N, D, n2);
if (testOut.amplitude>0.1)
drawRay(propogate(testOut, 5));
end
ylim([-10, 10])
xlim([-10, 10])
axis image
end
end
end
set(gcf, 'color', 'white')
axis image