

% clear memory
clear all;
clc;
% material defination
E = 1e-9; %input('Enter the value of Elastic Modulus : ');

v = 0.3; %input('Enter poission ratio : ');  %poission ratio

% matrix D
D=E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];

% load
P = 1000;
volfrac=0.6;
penal=4;
%Mesh generation
Lx=100; %input('Length of element at X axis : ');
Ly=30; %input('Length of element at Y axis : ');

nex = 100; %input('Number of element in X axis : ');
ney = 30; %input('Number of element in Y axis : ');
T=5;

numberElements=nex*ney;
[nodeCoordinates, elementNodes] = MeshMaker(Lx,Ly,nex,ney);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=2*numberNodes; 

% calculation of the system stiffness matrix
stiffness=calculate.stiffness(GDof,numberElements,...
 elementNodes,numberNodes,nodeCoordinates,D,T);


% solution
%displacement=zeros(2*(nex+1)*(ney+1),1);
force = zeros(2*(nex+1)*(ney+1),1);
force(2*(nex+1)*(ney+1)-nex+1,1)=-P;
force(2*(nex+1)*(ney+1)-nex-1,1)=-P;
force(2*(nex+1)*(ney+1)-nex-2,1)=-P;
%force((ney/2+1)(nex+1)-1,1)=3*P/2;
%force(2,1)=P/2;
%force(2*(nex+1)*(ney+1)-200,1)=P/2;
%force(1,1)=-P/4;
%force(2*(nex+1)*(ney+1)-1-200,1)=-P/4;
%stiffness4cal=stiffness;
bcdof=[2,4,2*(nex+1)-2,2*(nex+1)];
%bcdof=[1 2 7 8];
bcval=zeros(4,1);
[stiffness1,force1]=constraints(stiffness,force,bcdof,bcval);
displacements=stiffness1\force1;

displacements(isnan(displacements))=0;

force2=stiffness*displacements;

%stress calculation 

[Stress,Theta] = ElementalStress(nodeCoordinates,elementNodes,displacements,D);
vonmes=zeros(length(Stress),1);
for count=1:length(Stress)
vonmes(count,1)=sqrt((Stress(count,1))^2+(Stress(count,2))^2-(Stress(count,1))*(Stress(count,2))+3*abs((Stress(count,3))));

end

von=(reshape(vonmes,nex+1,ney+1))';

stressX=Stress(:,1);
stressXXX=(reshape(stressX',nex+1,ney+1))';

stressY=Stress(:,2);
stressYYY=(reshape(stressY',nex+1,ney+1))';

stressT=Stress(:,3);

stressTTT=(reshape(stressT',nex+1,ney+1))';

ThetaXYZ=reshape(Theta,nex+1,ney+1);
ThetaXY=ThetaXYZ';



 %Plot
 
figure;

imagesc(stressXXX);
title('stress along x axis');
axis off;
axis equal;
colormap(jet);
colorbar;
figure;

imagesc(stressYYY);
title('stress along y axis');
axis off;
axis equal;
colormap(jet);
colorbar;
figure;

imagesc(stressTTT);
title('shear stress');
axis off;
axis equal;
colormap(jet);
colorbar;


figure;
imagesc(von);
title('vonmess');
axis off;
axis equal;
colormap(jet);
colorbar;


