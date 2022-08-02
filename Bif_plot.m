function hs=Bif_plot(soln,p)
% Plot the bifurcation diagram in the (\hat a_1^{eff}, \hat a_2^{eff})
% parameter space and locate the bifurcation parameter values (returned by 
% handle hs) at top and rim of the colony.

% Locate node at top and rim
ind1=find(p.mesh.Nodes(2,:)==max(p.mesh.Nodes(2,:)));
ind2=find((abs(p.mesh.Nodes(2,:))< p.Eps)&(abs(p.mesh.Nodes(1,:)-p.R-p.U*p.dR)<p.Eps));
ind=[ind1,ind2];

% Read parameter and current solution values
u=soln(ind,1);
a1=p.A1;
a2=p.A2;
g1=[p.g1(1)*ones(length(ind1),1);(p.g1(1)+2*p.beta2)*ones(length(ind2),1)];
g2=[p.g2(1)*ones(length(ind1),1);(p.g2(1)+2*p.beta2)*ones(length(ind2),1)];
theta_x=p.theta_x;
theta_y=p.theta_y;
kp=p.kp;
km=p.km;

% Track C14/C4 for QS toggle and set the gate as 100% for NQS toggle
if p.QS
    C14=soln(ind,end-1)./(p.theta_g+soln(ind,end-1));
    C4=soln(ind,end)./(p.theta_h+soln(ind,end));
else
    C14=1;
    C4=1;
end

% Plot
symb=[repmat('o',1,length(ind1)) repmat('*',1,length(ind2))];
col=repmat('r',1,length(ind));
typ=[repmat('Top',length(ind1),1); repmat('Rim',length(ind2),1)];
plot(p.bif_a1,p.bif_a2,'k-','LineWidth',2);
hold on
hs=gscatter(a1./g1/theta_x.*C14,a2./g2/theta_y./(kp*u./(km+g2)+1).*C4,typ,col,symb,8);
xlim([0 6]);ylim([0 6]);
set(gca,'FontSize',20,'XTick',[0 3 6],'YTick',[0 3 6]);
axis square;
xlabel('$\hat a_1^{\textrm{eff}}$','Interpreter','latex');
ylabel('$\hat a_2^{\textrm{eff}}$','Interpreter','latex');
legend('Location','northeast');
hold off

end


