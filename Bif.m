function [bif_a1, bif_a2]=Bif(p)
syms x y a1 a2
param_max=20;
bif_a2=0:.01:param_max;
leak1=0;
leak2=0;
n1=2;n2=4;

bif_a1=zeros(2,length(bif_a2));
eq1=x-a1*leak1-a1*(1-leak1)/(1+y^n1);
eq2=y-a2*leak2-a2*(1-leak2)/(1+x^n2);
eq3=n1*n2*x^(n2-1)*(x-a1*leak1)^2*y^(n1-1)*(y-a2*leak2)^2-a1*(1-leak1)*a2*(1-leak2);
var=[x,y,a1];
eqs = [eq1, eq2, eq3];

for j=1:length(bif_a2)
    soln_range=[0,param_max;0 param_max;0.1 param_max;];
    eqs1=subs(eqs,a2,bif_a2(j));
    sol=vpasolve(eqs1,var,soln_range);
    if sol.a1
        bif_a1(:,j)=sort(sol.a1);
    else
        bif_a1(:,j)=Inf;
    end
end

end