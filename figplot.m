%% Plot
function figplot(model,soln,p)
u=soln(:,1);
x=soln(:,2);
y=soln(:,3);
row=2;colo=5;p1=[1 2 6 7];p2=[3 4 5];p3=[8 9 10]; SP=[p.po(1) p.po(2)+200 1500 500];
% row=2;colo=2;p1=[1 3];p2=2;p3=4;

fp=0;

if fp==1
   yt=soln(:,4); 
%    ybar=(yt-y)/2;
   row=3;
end
% rm=max(p.R_max,1);hm=max(p.H_max,0.2);

rm=3;hm=0.3;

hn=-1;
ind1=find(model.Mesh.Nodes(2,:)==max(model.Mesh.Nodes(2,:)));
list=[ind1,length(u)];
list=[];
dot_x=model.Mesh.Nodes(1,list);
dot_y=model.Mesh.Nodes(2,list);

cn=10;
YFP=ones(cn,3);YFP(:,1)=linspace(1,p.yellow(1),cn);YFP(:,2)=linspace(1,p.yellow(2),cn);YFP(:,3)=linspace(1,p.yellow(3),cn);
CFP=ones(cn,3);CFP(:,1)=linspace(1,p.blue(1),cn);CFP(:,2)=linspace(1,p.blue(2),cn);CFP(:,3)=linspace(1,p.blue(3),cn);

set(gcf,'Position',[0 200 1500 500]);

s(1)=subplot(row,colo,p1);
pdeplot(model,'XYData',u,'Contour','off','ColorMap',parula(cn*2));
text(dot_x,dot_y,num2str(u(list)));
u_max=max(p.u0,0.01);
caxis([0 u_max]);colorbar('Ticks',[0,u_max/2,u_max]);
axis equal;xlim([0 rm]);ylim([hn 0.5]);
%xlim([0.7 1]);ylim([0 0.03]);

title('aTc');%freezeColors;cbfreeze(colorbar);


s(2)=subplot(row,colo,p2);
pdeplot(model,'XYData',x,'Contour','off','ColorMap',YFP);
text(dot_x,dot_y,num2str(x(list)));
x_max=ceil(p.A1/p.g1);
caxis([0 x_max]);colorbar('Ticks',[0,x_max]);
axis equal;xlim([0 rm]);ylim([0 hm]);
title('LacI');%freezeColors;cbfreeze(colorbar);

s(3)=subplot(row,colo,p3);
pdeplot(model,'XYData',y,'Contour','off','ColorMap',CFP);
text(dot_x,dot_y,num2str(y(list)));
y_max=ceil(p.A2/p.g2);
caxis([0 y_max]);colorbar('Ticks',[0,y_max]);
axis equal;xlim([0 rm]);ylim([0 hm]);
title('TetR');%freezeColors;cbfreeze(colorbar);


if fp==1
    s(4)=subplot(row,colo,13:15);
    pdeplot(model,'XYData',yt,'Contour','off','ColorMap',CFP);
    text(dot_x,dot_y,num2str(yt(list)));
    caxis([0 5*p.theta_y]);
    axis equal;xlim([0 rm]);ylim([0 hm]);
    title('TetR^{Total}');%freezeColors;cbfreeze(colorbar);colormap(s(4),CFP);
end


colormap(s(1),parula(20));colormap(s(2),YFP);colormap(s(3),CFP);
set(s,'FontSize',18,'XTick', 0:1:rm, 'YTick', [hn:0.5:0 hm]);
set(gcf,'Position',SP,'color', 'w');
end