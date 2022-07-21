%% Plot
function figplot_tile(model,soln,p)
u=soln(:,1);
x=soln(:,2);
y=soln(:,3);

rm=max(p.R_max,1);
hm=max(p.H_max,0.2);
ind1=find(model.Mesh.Nodes(2,:)==max(model.Mesh.Nodes(2,:)));
list=[ind1,length(u)];

dot_x=model.Mesh.Nodes(1,list);
dot_y=model.Mesh.Nodes(2,list);

yellow=[255/255,211/255,0];blue=[75/255,134/255,180/255];cn=10;
YFP=ones(cn,3);YFP(:,1)=linspace(1,yellow(1),cn);YFP(:,2)=linspace(1,yellow(2),cn);YFP(:,3)=linspace(1,yellow(3),cn);
CFP=ones(cn,3);CFP(:,1)=linspace(1,blue(1),cn);CFP(:,2)=linspace(1,blue(2),cn);CFP(:,3)=linspace(1,blue(3),cn);

set(gcf,'Position',[0 200 1000 1000]);
tiledlayout(4,1);

s(1)=nexttile(1, [2 1]);
pdeplot(model,'XYData',u,'Contour','off','ColorMap',parula(cn*2));
text(dot_x,dot_y,num2str(u(list)));
caxis([0 p.u0]);
axis equal;xlim([0 rm]);ylim([-0.5 hm]);
%xlim([0.7 1]);ylim([0 0.03]);
title('aTc');%freezeColors;cbfreeze(colorbar);

s(2)=nexttile;
pdeplot(model,'XYData',x,'Contour','off','ColorMap',YFP);
text(dot_x,dot_y,num2str(x(list)));
caxis([0 3*p.theta_x]);
axis equal;xlim([0 rm]);ylim([0 hm]);
%xlim([0.7 1]);ylim([0 0.03]);
title('LacI');%freezeColors;cbfreeze(colorbar);

s(3)=nexttile;
pdeplot(model,'XYData',y,'Contour','off','ColorMap',CFP);
text(dot_x,dot_y,num2str(y(list)));
caxis([0 5*p.theta_y]);
axis equal;xlim([0 rm]);ylim([0 hm]);
%xlim([0.7 1]);ylim([0 0.03]);
title('TetR^{Free}');%freezeColors;cbfreeze(colorbar);


colormap(s(1),parula(20));colormap(s(2),YFP);colormap(s(3),CFP);
set(s,'FontSize',20,'TickLabelInterpreter','latex');
set(gcf,'Position',[0 200 1000 1000]);
end