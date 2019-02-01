    function [JacobianMatrix,Jacobianinverse]= Jacobian(nodeCoordinates,elementNodes,i,s,t)

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
   
    
    %defining jacobian
    J1=1/4*((-x1+x2+x3-x4)+(x1-x2+x3-x4)*t);
    J2=1/4*((-x1-x2+x3+x4)+(x1-x2+x3-x4)*s);
    J3=1/4*((-y1+y2+y3-y4)+(y1-y2+y3-y4)*t);
    J4=1/4*((-y1-y2+y3+y4)+(y1-y2+y3-y4)*s);
    
    JacobianMatrix=[J1 J3;J2 J4];                   

    Jacobianinverse=inv(JacobianMatrix);
    end 
