
function figplot_5(model,soln,p)
u=soln(:,1);
x=soln(:,2);
y=soln(:,3);
c14=soln(:,5);
c4=soln(:,6);

ind1=find(model.Mesh.Nodes(2,:)==max(model.Mesh.Nodes(2,:)));
ind=[ind1,length(u)];

dot_x=model.Mesh.Nodes(1,ind);
dot_y=model.Mesh.Nodes(2,ind);

% rm=max(p.R_max,1);hm=max(p.H_max,0.2);
rm=3;hm=0.3;

hn=-1;
cn=10;
YFP=ones(cn,3);YFP(:,1)=linspace(1,p.yellow(1),cn);YFP(:,2)=linspace(1,p.yellow(2),cn);YFP(:,3)=linspace(1,p.yellow(3),cn);
CFP=ones(cn,3);CFP(:,1)=linspace(1,p.blue(1),cn);CFP(:,2)=linspace(1,p.blue(2),cn);CFP(:,3)=linspace(1,p.blue(3),cn);

Nr=4;Nc=6;

s(1)=subplot(Nr,Nc,[1,2,7,8]);
pdeplot(model,'XYData',u,'Contour','off','ColorMap',parula(cn*2));
text(dot_x,dot_y,num2str(u(ind),'%.2f'));
u_max=max(p.u0,0.01);
caxis([0 u_max]);colorbar('Ticks',[0,u_max/2,u_max]);
axis equal;xlim([0 rm]);ylim([hn hm]);
title('aTc');%freezeColors;cbfreeze(colorbar);


s(2)=subplot(Nr,Nc,13:18);
pdeplot(model,'XYData',x,'Contour','off','ColorMap',YFP);
text(dot_x,dot_y,num2str(x(ind),'%.2f'));
axis equal;xlim([0 rm]);ylim([0 hm]);
title('LacI');%freezeColors;cbfreeze(colorbar);
x_max=max(ceil(p.A1/p.g1),1);
caxis([0 x_max]);colorbar('Ticks',[0,x_max/2,x_max]);

s(3)=subplot(Nr,Nc,19:24);
pdeplot(model,'XYData',y,'Contour','off','ColorMap',CFP);
text(dot_x,dot_y,num2str(y(ind),'%.2f'));
axis equal;xlim([0 rm]);ylim([0 hm]);
title('TetR');%freezeColors;cbfreeze(colorbar);
y_max=max(ceil(p.A2/p.g2),1);
caxis([0 y_max]);colorbar('Ticks',[0,y_max/2,y_max]);

s(4)=subplot(Nr,Nc,[3,4,9,10]);
pdeplot(model,'XYData',c14,'Contour','off','ColorMap',YFP);
c14_max=max(ceil(p.theta_g*10),1);
caxis([0 c14_max]);colorbar('Ticks',[0,c14_max/2,c14_max]);
text(dot_x,dot_y,num2str(c14(ind),'%.2f'));
axis equal;xlim([0 rm]);ylim([hn hm]);
title('C14');%freezeColors;cbfreeze(colorbar);

s(5)=subplot(Nr,Nc,[5,6,11,12]);
pdeplot(model,'XYData',c4,'Contour','off','ColorMap',CFP);
c4_max=max(ceil(p.theta_h*10),1);
caxis([0 c4_max]);colorbar('Ticks',[0,c4_max/2,c4_max]);
text(dot_x,dot_y,num2str(c4(ind),'%.2f'));
axis equal;xlim([0 rm]);ylim([hn hm]);
title('C4');%freezeColors;cbfreeze(colorbar);

colormap(s(1),parula(cn*2));
colormap(s(2),YFP);colormap(s(3),CFP);
colormap(s(4),YFP);colormap(s(5),CFP);
set(s,'FontSize',18,...
      'XTick', 0:1:rm, 'YTick', [hn:0.5:0 hm]);
set(s(4:5),'XTick',{});
set(gcf,'Position',[p.po(1) p.po(2)+200 min(1200,p.po(3)) min(800,p.po(4))],'color', 'w');
end