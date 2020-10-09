# calculation of error in displacement
function [u] =Err(s,f,u)
e=s*u-f;
error =abs(mean(e));
while error > 0.01
ue=s\e;
u1=u
f1=f-e;
u1=s\f1;
e=f-s*u1;
error=abs(mean(e));
end
u=u1;
end

