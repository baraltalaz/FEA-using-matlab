function [stiffness,force]=constraints(stiffness,force,constrainsdof,constrainsval)
 n=length(constrainsdof);
 sdof=length(stiffness);

 for i=1:n
    c=constrainsdof(i);
    for j=1:sdof
       stiffness(c,j)=0;
    end

    stiffness(c,c)=1;
    force(c)=constrainsval(i);
 end
