function hs=Bif_plot(soln,p)

ind1=find(p.mesh.Nodes(2,:)==max(p.mesh.Nodes(2,:)));
ind2=find((abs(p.mesh.Nodes(2,:))< p.Eps)&(abs(p.mesh.Nodes(1,:)-p.R-p.U*p.dR)<p.Eps));

% ind2=find((abs(p.mesh.Nodes(2,:))< 1e-5)&(p.mesh.Nodes(1,:)<p.R+(p.U-1)*p.dR));
% ind=[ind1,length(soln(:,1))];
ind=[ind1,ind2];
if length(ind)<2
    stop
end

u=soln(ind,1);

a1=p.A1;
a2=p.A2;
% g1=[p.g1;p.g1+p.beta2;p.g1+2*p.beta2];
g1=[p.g1(1)*ones(length(ind1),1);(p.g1(1)+2*p.beta2)*ones(length(ind2),1)];
g2=[p.g2(1)*ones(length(ind1),1);(p.g2(1)+2*p.beta2)*ones(length(ind2),1)];%%%%%%%%%%%%%%%%%%%%%%
theta_x=p.theta_x;
theta_y=p.theta_y;
kp=p.kp;
km=p.km;

if p.QS
    C14=soln(ind,end-1)./(p.theta_g+soln(ind,end-1));
    C4=soln(ind,end)./(p.theta_h+soln(ind,end));
else
    C14=1;
    C4=1;
end
symb=[repmat('o',1,length(ind1)) repmat('*',1,length(ind2))];
col=repmat('r',1,length(ind));
typ=[repmat('Top',length(ind1),1); repmat('Rim',length(ind2),1)];
plot(p.bif_a1,p.bif_a2,'k-','LineWidth',2);
hold on
% hs=gscatter(a1./g1/theta_x.*C14,a2./g2/theta_y./(kp*u./(km+g2)+1).*C4,{'Top';'AGL';'Tip'},'rrr','ox*',8);
hs=gscatter(a1./g1/theta_x.*C14,a2./g2/theta_y./(kp*u./(km+g2)+1).*C4,typ,col,symb,8);
% set(hs(:),'Color',col(i,:));
% xlim([0 10]);ylim([0 10]);
    xlim([0 6]);ylim([0 6]);
    set(gca,'FontSize',20,'XTick',[0 3 6],'YTick',[0 3 6]);

% if p.QS
%     xlim([0 6]);ylim([0 6]);
%     set(gca,'FontSize',20,'XTick',[0 3 6],'YTick',[0 3 6]);
% else
%     xlim([0 3]);ylim([0 6]);
%     set(gca,'FontSize',20,'XTick',[0 3],'YTick',[0 3 6]);
% end
% xlim([1 3]);ylim([1 3]);

set(gca,'FontSize',20,'XTick',[0 3 6],'YTick',[0 3 6]);
axis square;
xlabel('$\hat a_1^{\textrm{eff}}$','Interpreter','latex');
ylabel('$\hat a_2^{\textrm{eff}}$','Interpreter','latex');
legend('Location','northeast');
% set(gca,'linewidth',2);
%     f=get(gca,'Children');
%     legend([f(3),f(5),f(1)]);
hold off

end


