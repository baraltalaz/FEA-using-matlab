classdef calculate
   methods (Static)
     % stiffness calculation...............
        function K = stiffness(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,D,thickness)
            
         K=sparse(zeros(GDof));   
         [gaussWeights,gaussLocations]=calculate.GaussQuadrature('Full');
         
         for i=1:numberElements                           
  ElementNo=elementNodes(i,:); 
  elementDof=[2*ElementNo(1)-1 , 2*ElementNo(1)   2*ElementNo(2)-1   2*ElementNo(2)   2*ElementNo(3)-1   2*ElementNo(3)  2*ElementNo(4)-1 2*ElementNo(4)]; 
  
  ndof=length(ElementNo);
  
  % Gauss point
  B =zeros(8);
  for j=1:size(gaussWeights,1)                      
    GaussPoint=gaussLocations(j,:);                                                     
    s=GaussPoint(1);
    t=GaussPoint(2);
    
% shape functions and natural derivatives
   % [shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta);
shapeFunction=1/4*[ (1-s)*(1-t);(1+s)*(1-t);(1+s)*(1+t);(1-s)*(1+t)];
ND = 1/4*[-(1-t), -(1-s);1-t, -(1+s);1+t, 1+s;-(1+t), 1-s];
% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
ElementNodes=elementNodes(i,:);
[Jacob,Jacobianinverse]=Jacobian(nodeCoordinates,ElementNodes,i,s,t);

%  B matrix
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
    
    
    B=det(Jacobianinverse)*b;
    
    
% stiffness matrix
    K(elementDof,elementDof)=...
        K(elementDof,elementDof)+...
        B'*D*thickness*B*gaussWeights(j)*det(Jacob);    

    
    
  end  
end    

         
         
         
        end
      % Gauss weight calling function ...............  
       function [weights,locations] = GaussQuadrature(option)
            
         switch option
        case 'Full'
    
        locations=...
          [ -0.577350269189626 -0.577350269189626;
             0.577350269189626 -0.577350269189626;
             0.577350269189626  0.577350269189626;
            -0.577350269189626  0.577350269189626];
        weights=[ 1;1;1;1]; 
    
        case 'reduced'
        
        locations=[0 0];
        weights=[4];
    end
        end
    
    
    
    end
end


