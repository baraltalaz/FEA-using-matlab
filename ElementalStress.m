function [Stress,Theta] = ElementalStress(nodeCoordinates,elementNodes,U,D)

Stress=zeros(length(nodeCoordinates),3);
Theta= zeros(length(nodeCoordinates),1);
for i= 1:size(elementNodes(:,1))
    %{
    x1 = nodeCoordinates(elementNodes(i,1),1);
    y1 = nodeCoordinates(elementNodes(i,1),2);
    x2 = nodeCoordinates(elementNodes(i,2),1);
    y2 = nodeCoordinates(elementNodes(i,2),2);
    x3 = nodeCoordinates(elementNodes(i,3),1);
    y3 = nodeCoordinates(elementNodes(i,3),2);
    x4 = nodeCoordinates(elementNodes(i,4),1);
    y4 = nodeCoordinates(elementNodes(i,4),2);
    %}
    
    a=elementNodes(1,1);
    b=elementNodes(1,2);
    c=elementNodes(1,3);
    d=elementNodes(1,4);
    x1 = nodeCoordinates(a,1);
    y1 = nodeCoordinates(a,2);
    x2 = nodeCoordinates(b,1);
    y2 = nodeCoordinates(b,2);
    x3 = nodeCoordinates(c,1);
    y3 = nodeCoordinates(c,2);
    x4 = nodeCoordinates(d,1);
    y4 = nodeCoordinates(d,2);
    u=[U(2*a-1) U(2*a) U(2*b-1) U(2*b) U(2*c-1) U(2*c) U(2*d-1) U(2*d)];
    UU=u';
    [gaussWeights,gaussLocations]=calculate.GaussQuadrature('Full');
    
    X=[x1 x2 x3 x4];
    Y=[y1 y2 y3 y4];
        
   B=zeros(3,8);
   
     for j=1:size(gaussWeights,1)                      
    GaussPoint=gaussLocations(j,:);                                                     
    s=GaussPoint(1);
    t=GaussPoint(2);
    
      shapeFunction=1/4*[ (1-s)*(1-t);(1+s)*(1-t);(1+s)*(1+t);(1-s)*(1+t)];
     ND = 1/4*[-(1-t), -(1-s);1-t, -(1+s);1+t, 1+s;-(1+t), 1-s];
   
     ElementNodes=elementNodes(i,:);
     
     [Jacob,Jacobianinverse]=Jacobian(nodeCoordinates,ElementNodes,i,s,t);
%{
 B matrix
    b=zeros(3,8);
    b(1,1) = ND(1,1); b(1,3) = ND(2,1); b(1,5) = ND(3,1); b(1,7) = ND(4,1);
    b(2,2) = ND(1,2);b(2,4) = ND(2,2); b(2,6) = ND(3,2); b(2,8) = ND(4,2);
     b(3,1) = ND(1,2);
    b(3,2) = ND(1,1);
    b(3,3) = ND(2,2);
    b(3,4) = ND(2,1);
    b(3,5) = ND(3,2);
    b(3,6) = ND(3,1);
    b(3,7) = ND(4,2);
    b(3,8) = ND(4,1);
    %}
     x=X(j);
     y=Y(j);
    b=[y-y4 0 y3-y 0 y-y2 0 y1-y 0 ;0 x-x2 0 x1-x 0 x-x4 0 x3-x ; x-x2 y-y4 x1-x y3-y x-x4 y-y2 x3-x y1-y ];        
    
   B=det(Jacobianinverse)*b;
    s=D*B*UU;
    S=s';
    posi=ElementNodes(j);
    Stress(posi,:)=Stress(posi,:)+S;
   Theta(posi,1)=(1/2)*(180/pi)*(atan((2*Stress(i,3))/(Stress(i,1)-Stress(i,2))));
end   
 
  end
end

