function circleplot(p,soln)
% Plot the subplots of the top-down views of the colony with LacI profile, 
% TetR profile and the merged profiles.

loc_c=round(p.mesh.Nodes,5);

% Record colony boundary (r-z space) 
xv=[0, 0, p.R+p.U*p.dR];
yv=[0, max(loc_c(2,:)), 0];

% Generating a finer colony mesh
dr=p.dmesh/10;
r = 0:dr:(p.R+p.U*p.dR);
h = 0:dr:max(loc_c(2,:)); 
[R,H] = meshgrid(r,h);

% Projecting solution to finer colony mesh
in = inpolygon(R,H,xv,yv);

Z_LacI=0*R;Z_TetR=Z_LacI;

F1 = scatteredInterpolant(loc_c(1,:)',loc_c(2,:)',soln(:,2));
Z_LacI(in) = F1(R(in),H(in));

F2 = scatteredInterpolant(loc_c(1,:)',loc_c(2,:)',soln(:,3));
Z_TetR(in) = F2(R(in),H(in));


% Integrating the solution in the r-direction
Int_Z_LacI=0*r(1:end-1);Int_Z_TetR=Int_Z_LacI;
for i=1:(length(r)-1)
    Int_Z_LacI(i)=trapz(r(i:i+1),trapz(h,Z_LacI(:,i:i+1),1));
    Int_Z_TetR(i)=trapz(r(i:i+1),trapz(h,Z_TetR(:,i:i+1),1));
end

% Define colormap
cmap_yfp = (0:0.01:1)'*p.yellow_bright; 
cmap_cfp = (0:0.01:1)'*p.blue_bright; 

% Define range of colorbar
x_max_circ=p.x_max*dr*p.H_max/2;
y_max_circ=p.y_max*dr*p.H_max/2;

% Genrating mesh and solution for the top-down (circle) view
model = createpde;
geometryFromEdges(model,decsg([1;0;0;(p.R+p.U*p.dR)]));
m1 = generateMesh(model,'Hmax',0.01);
x=m1.Nodes(1,:);
y=m1.Nodes(2,:);
r_mesh=sqrt(x.^2+y.^2);
sol_laci=interp1(r(1:end-1),Int_Z_LacI,r_mesh);
sol_TetR=interp1(r(1:end-1),Int_Z_TetR,r_mesh);
range_plot=[-p.rm p.rm]; 

% Adding the top-down view suplots: LacI, TetR, Merged
s(1)=subplot(4,4,[9 13]);
pdeplot(model,'XYData',sol_laci,'Contour','off');
axis equal;
set(gca,'color','k','FontSize',20);
caxis([0 x_max_circ]);
colorbar off;
colorbar('south','Ticks',[0,x_max_circ],'TickLabels',{0,1},'Position',[0.1300    0.1    0.1569    0.0160]);
title('LacI');

s(2)=subplot(4,4,[10 14]);
pdeplot(model,'XYData',sol_TetR,'Contour','off');
axis equal;
set(gca,'color','k','FontSize',20);
caxis([0 y_max_circ]);
colorbar off;
colorbar('south','Ticks',[0,y_max_circ],'TickLabels',{0,1},'Position',[0.3362    0.1    0.1569    0.0160]);
title('TetR');

colormap(s(1),cmap_yfp);
colormap(s(2),cmap_cfp);

% Top-down Merged
s(3)=subplot(4,4,[11 15]);
t=model.Mesh.Elements;
P=model.Mesh.Nodes;
it1=t(1,:);
it2=t(2,:);
it3=t(3,:);
X=[P(1,it1); P(1,it2); P(1,it3)];
Y=[P(2,it1); P(2,it2); P(2,it3)];
Cx=[sol_laci(it1);sol_laci(it2);sol_laci(it3)];      
Cy=[sol_TetR(it1);sol_TetR(it2);sol_TetR(it3)];

% Define 2D gradient
C=zeros([size(X),3]);
C(:,:,1)=p.yellow_bright(1)*Cx/x_max_circ+p.blue_bright(1)*Cy/y_max_circ;
C(:,:,2)=p.yellow_bright(2)*Cx/x_max_circ+p.blue_bright(2)*Cy/y_max_circ;
C(:,:,3)=p.yellow_bright(3)*Cx/x_max_circ+p.blue_bright(3)*Cy/y_max_circ;
patch(X,Y,0.*X,C,'Edgecolor','none');
title('Merged');
axis equal;
set(gca,'color','k','FontSize',20);

set(s,'xlim',range_plot,'ylim',range_plot,'xtick',[-p.rm 0 p.rm],'ytick',[-p.rm 0 p.rm]);
end
