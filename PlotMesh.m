# matlab codes to plot mesh
function PlotMesh(coordinates,nodes)
nel = length(nodes) ;                  % number of elements
nnode = length(coordinates) ;          % total number of nodes in system
nnel = size(nodes,2);                % number of nodes per element
% 
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;

for iel=1:nel   
     for i=1:nnel
     nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
     X(i,iel)=coordinates(nd(i),1);    % extract x value of the node
     Y(i,iel)=coordinates(nd(i),2);    % extract y value of the node
     end
end
    
% Plotting the FEM mesh, diaplay Node numbers and Element numbers
     f1 = figure ;
     set(f1,'name','Mesh','numbertitle','off') ;
     plot(X,Y,'k')
     %fill(X,Y,'w')
     colormap(jet);
     title('Finite Element Mesh') ;
     axis off ;
     k = nodes(:,1:end);
     nd = k' ;
    for i = 1:nel
     text(X(:,i),Y(:,i),int2str(nd(:,i)),'fontsize',8,'color','k');
     text(sum(X(:,i))/4,sum(Y(:,i))/4,int2str(i),'fontsize',10,'color','r') ;
    end        
