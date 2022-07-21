clear
% close all;
Tmax=15;
global theta_x theta_y A1 A2 g1 g2 kp km A3 A4 theta_g theta_h g_C14 g_C4 p
p.dT=10;
p.T=15;
gamma_1=log(2)/7*p.dT;
dil=log(2)/p.T*log(4)/log(2);
A1=100*p.dT;A2=200*p.dT;
g1=gamma_1+0*dil;
g2=g1;
theta_x=500;theta_y=500;
p.l1=0.2;
kp=0.06*p.dT/200;km=kp/10;

g_C14=log(2)/(24*6)+dil;
g_C4=g_C14;
A3=.40;A4=.80;
p.l3=0.75; p.l4=0.1;

theta_g=0;theta_h=0;

u_range=10.^(0:0.1:3);
% u_range=23.5;
Yall=zeros(length(u_range),2);
for i = 1:length(u_range)
    u=u_range(i);
    [t,Y]=ode45(@(t,X) ode_QS(t,X,u),[0,Tmax],[0,0,0,2e5,2e5]);
    Y_all(i,:)=Y(end,1:2);
end

figure,
semilogx(u_range*0.4255,Y_all(:,1),"Color",[0.9290, 0.6940, 0.1250]);
hold on
semilogx(u_range*0.4255,Y_all(:,2),"Color",[0, 0.4470, 0.7410]);
hold off
xlim([1 100]);


fig=figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

% yyaxis left
h1=plot(10*t,[Y(:,1),Y(:,2)],'-','LineWidth',3);
set(h1, {'color'}, {[0.9290, 0.6940, 0.1250]; [0, 0.4470, 0.7410];});
ylabel('$x(t),\; y(t)$','Interpreter','latex')
% ylim([0 max([Y(:,1);Y(:,2)])])

% yyaxis right
% h2=plot(10*t,[Y(:,4),Y(:,5)],'--','LineWidth',3);
% set(h2, {'color'}, {[0.9290, 0.6940, 0.1250]; [0, 0.4470, 0.7410];});
% ylabel('$g(t),\; h(t)$','Interpreter','latex')

set(gca,'FontSize',20,'TickLabelInterpreter','latex')
legend('LacI','TetR','C14','C4','Interpreter','latex')
xlabel('Time $t$ (min)','Interpreter','latex')
title('$a_1=2; a_2=5; u_0=2$','Interpreter','latex');

function rhs=ode_QS(t,X,u)
global theta_x theta_y A1 A2 g1 g2 A3 A4 theta_g theta_h kp km g_C14 g_C4 p

rhs=X;
x=X(1);
y=X(2);
yt=X(3);
g=X(4);
h=X(5);
alpha=1/4;
rhs(1)= A1*p.l1+A1*(1-p.l1)/(1+(y/theta_y)^2).*(g/(theta_g+g))-g1*x;
rhs(2)= A2/(1+(x/theta_x)^4).*h/(theta_h+h)-kp*u.*y+km*(yt-y)-g2*y;
rhs(3)= A2/(1+(x/theta_x)^4).*h/(theta_h+h)-g2*yt;
rhs(4)= A3*p.l3+A3*(1-p.l3)/(1+(y/theta_y)^2)-g_C14*g;
rhs(5)= A4*p.l4+A4*(1-p.l4)/(1+(x/theta_x)^4)-g_C4*h;
end