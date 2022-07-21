function ic=Initialization(mesh_old,mesh,soln,p)

dH_m=p.dH_mesh;
nodes_old=round(mesh_old.Nodes,p.n_eps);
R1=p.R+(p.U-p.dU)*p.dR;

% 1. Push all nodes in the Top domain up by dH
% N_top=findNodes(mesh_old,'region','Face',3);
N_top=find(nodes_old(2,:)>=(p.H_agl-p.Eps));
nodes_top=nodes_old(:,N_top);
nodes_top(2,:)=nodes_top(2,:)+dH_m; 

% 2. Stretch the shape of the AGL+Rim
% N_agl=findNodes(mesh_old,'region','Face',2);
N_agl=find(nodes_old(1,:)<=(R1+p.Eps) & abs(nodes_old(2,:)-p.H_agl/2)<=(p.H_agl/2+p.Eps));
nodes_agl=nodes_old(:,N_agl);

% 2.1 Stretch the z-coordinates of AGL+Rim by (1+dH_m/p.dH)
ind_l=nodes_old(2,N_agl)>p.Eps;
nodes_agl(2,ind_l)=nodes_agl(2,ind_l)*(1+p.dU);

% 2.2 Stretch the r-coordinates of the Rim domain by (1+dR_m/p.dR)
R2=p.R+(p.U-1-p.dU)*p.dR;
ind_r=nodes_old(1,N_agl)>(R2-p.Eps);
nodes_agl(1,ind_r)=R2+(nodes_agl(1,ind_r)-R2)*(1+p.dU);

nodes_old(:,N_top)=nodes_top;
nodes_old(:,N_agl)=nodes_agl;
nodes_old=round(nodes_old,p.n_eps);

% % 3. Projecting the solution to a finer uniform mesh
% p.dmesh=p.dmesh/10;
% [~, mesh_finer]=MeshGenerator(p,1);
% nodes_finer = mesh_finer.Nodes;
% N1=findNodes(mesh_finer,'region','Face',[2 3]);
% N2=findNodes(mesh_finer,'region','Edge',[4:6 9 10]);
% N0=N1(~ismember(N1,N2));
% 
% nodes_finer(:,N0)=[];
% 
% H=max(nodes_old(2,:));
% R=p.R+p.U*p.dR;
% 
% xv=[0, 0, R];
% yv=[0, H, 0];
% 
% r = 0:p.dmesh:R;
% h = 0:p.dmesh:H;
% 
% [Rq,Hq] = meshgrid(r,h);
% [in, on] = inpolygon(Rq,Hq,xv,yv);
% RHq=[Rq(in&~on),Hq(in&~on)]';
% 
% nodes_finer=round([nodes_finer,RHq],p.n_eps);
% 
% M=size(nodes_finer,2);
% ic_finer=zeros(M,p.N);
% 
% for i=1:p.N
%     F_finer =scatteredInterpolant(nodes_old(1,:)',nodes_old(2,:)',soln(:,i));
%     if ismember(i,2:4)
%         F_finer.Method = 'nearest';
%     else
%         F_finer.Method = 'natural';
%     end
%     ic_finer(:,i) = F_finer(nodes_finer(1,:),nodes_finer(2,:));
% end

nodes_finer=nodes_old;
ic_finer=soln;
% 4. Coarsening the finer interpolated solution from previous step to the new mesh
nodes_new=round(mesh.Nodes,p.n_eps);
M0=size(nodes_new,2);

ic=zeros(M0,p.N);

for i=1:p.N
    F =scatteredInterpolant(nodes_finer(1,:)',nodes_finer(2,:)',ic_finer(:,i));
%     F =scatteredInterpolant(nodes_finer(1,:)',nodes_finer(2,:)',ic_finer(:,i));
    if ismember(i,2:4)
        F.Method = 'linear';
        F.ExtrapolationMethod = 'none';
    else
        F.Method = 'natural';
        F.ExtrapolationMethod = 'none';
    end
    ic(:,i) = F(nodes_new(1,:),nodes_new(2,:));
end

top=(nodes_new(2,:)>-p.Eps)&(nodes_new(1,:)<(p.R+p.U*p.dR));
ic(~top,2:4)=0;

ID=isnan(ic);
if sum(ID(:))
    stop
end

% figure
% pdeplot(mesh)
% hold on
% h1=text(nodes_new(1,:),nodes_new(2,:),num2str(ic(:,2)));
% axis equal; xlim([p.R+(p.U-1)*p.dR p.R+p.U*p.dR]); ylim([0 p.H_agl]);
% hold off
% set(gcf,'Position',[p.po(1) p.po(2) p.po(3) 400], 'color', 'w','Resize','off');
% 
% figure
% pdeplot(mesh_old)
% hold on
% h2=text(mesh_old.Nodes(1,:),mesh_old.Nodes(2,:),num2str(soln(:,2)));
% axis equal; xlim([R2 R1]); ylim([0 p.H_agl]);
% hold off
% set(gcf,'Position',[p.po(1) p.po(2)+400 p.po(3) 400], 'color', 'w','Resize','off');

end


%% Test
% pdeplot(mesh_old);
% hold on
% plot(nodes_old(1,:),nodes_old(2,:),'*')
% plot(mesh_old.Nodes(1,N_top),mesh_old.Nodes(2,N_top)+0.001,'o')
% plot(mesh_old.Nodes(1,N_agl(ind_l)),mesh_old.Nodes(2,N_agl(ind_l))*1.1,'^')
% plot(R2+(mesh_old.Nodes(1,N_agl(ind_r))-R2)*1.1,mesh_old.Nodes(2,N_agl(ind_r))*1.1,'s')


% tic
% model = createpde(p.N);
% geometryFromMesh(model,nodes_old,mesh_old.Elements,ones(1,size(mesh_old.Elements,2)));
% soln_m=createPDEResults(model,reshape(soln,[],1));
% ic=interpolateSolution(soln_m,nodes_new,1:(p.N));
% N_agar=findNodes(mesh,'region','Face',1);
% ic(N_agar,2:4)=0;
% toc


% figure,pdeplot(mesh_old);h1=text(mesh_old.Nodes(1,N0),mesh_old.Nodes(2,N0),num2str(soln(N0,i)));
% figure, pdeplot(mesh); h1=text(mesh.Nodes(1,N1),mesh.Nodes(2,N1),num2str(ic(N1,3),1)); xlim([0.48 0.52]); ylim([0.06 0.08]);
% set(gcf, 'Position', get(0, 'Screensize'));



% figure,
% i=3;
% % EdgeID=[4,5];
% 
% N_edge0=findNodes(mesh_old,'region','Edge',EdgeID);
% plot(nodes_old(1,N_edge0),soln(N_edge0,i),'-','LineWidth',3,'MarkerSize',10);
% N_edge1=findNodes(mesh,'region','Edge',EdgeID);
% hold on
% plot(mesh.Nodes(1,N_edge1),ic(N_edge1,i),'o','LineWidth',3,'MarkerSize',10);
% legend({'Before', 'After'},'Location','northwest');
% set(gca,'FontSize',24);
% hold off
