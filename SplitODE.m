function rhs=SplitODE(~,X,p,indx_i)

rhs = zeros(3,length(indx_i));
X = reshape(X,[],length(indx_i));

ic=p.ic';

u=ic(1,:);
x =X(1,:);
y =X(2,:);
yt=X(3,:);
ybar=yt-y;

if p.QS==0
    prod_x=p.A1./(1+(y/p.theta_y).^2);
    prod_y=p.A2./(1+(x/p.theta_x).^4);
elseif p.QS==1
    g=ic(p.N-1,indx_i);
    h=ic(p.N,indx_i);
    prod_x=p.A1*p.leak1+p.A1.*(1-p.leak1)./(1+(y/p.theta_y).^2).*(g./(p.theta_g+g));
    prod_y=p.A2./(1+(x/p.theta_x).^4).*h./(p.theta_h+h);
 end
sqst=-p.kp*u(indx_i).*y+p.km*ybar;
rhs(1,:) =prod_x-(p.g1+p.dilution(indx_i)).*x;
rhs(2,:) =prod_y-(p.g2+p.dilution(indx_i)).*y+sqst;
rhs(3,:) =prod_y-(p.g2+p.dilution(indx_i)).*yt;

rhs = rhs(:);

end