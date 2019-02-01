function [coordinates,nodes] = MeshMaker(L,B,Nx,Ny)
x=0; y=0;
xinc=L/Nx;
yinc=B/Ny;
coordinates=[];
 for j=1:Ny+1
     for i=1:Nx+1
        coordinates=[coordinates;x y];
         x=x+xinc;
      end
     x=0;
     y=y+yinc;
 end
 
yth=0;
  nel=0;
  count=0;
for k=1:Ny
for l=1:Nx
    nel=yth*(Nx+1)+l;
    count=count+1;
        nodes(count,:)=[nel nel+1  nel+Nx+2 nel+Nx+1];
        
        
end
yth=yth+1;
end

