function ic=Initialization(mesh_old,mesh,soln,p)
% Initialize the new nodes born from AGL+Rim

dH_m=p.dH_mesh;
nodes_old=round(mesh_old.Nodes,p.n_eps);
R1=p.R+(p.U-p.dU)*p.dR;

% 1. Push all nodes in the Top domain up by dH
N_top=find(nodes_old(2,:)>=(p.H_agl-p.Eps));
nodes_top=nodes_old(:,N_top);
nodes_top(2,:)=nodes_top(2,:)+dH_m; 

% 2. Stretch the shape of the AGL+Rim
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


% 3. Interpolated solution from previous step to the new mesh
nodes_new=round(mesh.Nodes,p.n_eps);
M0=size(nodes_new,2);

ic=zeros(M0,p.N);

for i=1:p.N
    F =scatteredInterpolant(nodes_old(1,:)',nodes_old(2,:)',soln(:,i));
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

end
