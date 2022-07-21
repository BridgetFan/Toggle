clear
% close all;
Tmax=50;
global theta_x theta_y A1 A2 g1 g2 kp km 

dT=10;
A1=75*dT;
A2=150*dT;
g1=.1*dT;
g2=g1;
theta_x=500;
theta_y=500;
kp=0.06*dT;
km=kp/10;
alpha=0.8;


u=30;
[t,Y]=ode45(@(t,X) ode_NQS(t,X,u),[0,Tmax],[0,0,0]);
figure
h1=plot(10*t,[Y(:,1),Y(:,2)],'LineWidth',3);
set(h1, {'color'}, {[0.9290, 0.6940, 0.1250]; [0, 0.4470, 0.7410];});
set(gca,'FontSize',20,'TickLabelInterpreter','latex')
xlabel('Time $t$ (min)','Interpreter','latex')
ylabel('$x(t),\; y(t)$','Interpreter','latex')
legend('LacI','TetR','Interpreter','latex')
ylim([0 max([Y(:,1);Y(:,2)])])
tit=['$\hat a_1=' num2str(A1/g1/theta_x) ', \hat a_2=' num2str(A2/g2/theta_y) ', u_0=' num2str(u) '$'];
title(tit,'Interpreter','latex');


function rhs=ode_NQS(t,X,u)
global theta_x theta_y A1 A2 g1 g2 kp km 
rhs=X;
x=X(1);
y=X(2);
yt=X(3);

rhs(1)= A1/(1+(y/theta_y)^2)-g1*x;
rhs(2)= A2/(1+(x/theta_x)^4)-kp*u.*y+km*(yt-y)-g2*y;
rhs(3)= A2/(1+(x/theta_x)^4)-g2*yt;
end