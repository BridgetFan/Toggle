function [model, mesh]=MeshGenerator(p)

model = createpde(p.N);

R1=round(p.R+p.U*p.dR,p.n_eps);
H1=round(p.H+p.U*p.dH,p.n_eps);

R2=round(R1-p.dR,p.n_eps);
H2=round(p.H_agl,p.n_eps);

% dz=p.dH*p.dt_mesh/p.T;

xc=round(p.dR:p.dR:R2,p.n_eps);
yc=round((H1-p.dH):-p.dH:H2,p.n_eps);

if mod(p.U,1)
    xc=[xc R2];
    yc=[yc H2];
end

% xc=[xc,R2-p.dH:-p.dH:p.dH];
% yc=[yc,H2+p.dH/10:p.dH/10:H1-p.dH/10];
% Z=ones(size(xc));
% x_list=[0*Z;0*Z;xc];
% y_list=[H1*Z;yc;yc];
% Corrd=[x_list(:),y_list(:)];
%
% Top = [2  3  0   R2  0 ...
%              H1  H2  H2]';


x_list=round([0 0:p.dH:R2]',p.n_eps);
X=repmat(x_list,1,length(xc));
X1=repmat(xc,size(X,1),1);
ind=X>X1;
X(ind)=nan;
Y=repmat(yc,size(X,1),1);
Y(ind)=nan;
Y(1,:)=H1;
Corrd=[X(~isnan(X(:))) Y(~isnan(Y(:)))];

Top = [2 size(Corrd,1) Corrd(:)']';

% AGL = [2  4  0   0  R1  R2 ...
%              H2  0  0   H2]';
AGL = [2  length(0:p.dH:R1)+2  0   0:p.dH:R1  R2 ...
    H2  zeros(size(0:p.dH:R1))   H2]';

Agar = [3  4  0  p.Ra  p.Ra   0 ...
    0  0    -p.Ha  -p.Ha]';

L_max=max([length(Top),length(AGL),length(Agar)]);
Top = [Top;zeros(L_max - length(Top),1)];
AGL = [AGL;zeros(L_max - length(AGL),1)];
Agar = [Agar;zeros(L_max - length(Agar),1)];

% gm = [Top,AGL,Agar];
% sf = 'Top+AGL+Agar';
% ns = char('Top','AGL','Agar');

gm = round([Agar,AGL,Top],p.n_eps);
sf = 'Agar+AGL+Top';
ns = char('Agar','AGL','Top');

ns = ns';
g = decsg(gm,sf,ns);
pg=geometryFromEdges(model,g);

% pdegplot(pg, 'VertexLabels', 'off', 'EdgeLabels', 'off','FaceLabels','on','FaceAlpha',0.5);

% edgeIDs = faceEdges(pg,1,'external');
% mesh = generateMesh(model,'Hmax',p.dmesh,'Hedge',{edgeIDs,0.2},'Hface',{1,0.2});


% pdegplot(model,'EdgeLabels','on','FaceLabels','on');
% axis equal
% xlim([-1.1,1.1])
% m1 = generateMesh(model);

% mesh = generateMesh(model,'Hedge',{[3:6,9,10],p.dmesh});%,'GeometricOrder','linear');%%%%%%%%%%%%%%
% mesh = generateMesh(model,'Hface',{2:3,p.dmesh});%,'GeometricOrder','linear');%%%%%%%%%%%%%%

F_t=ceil(round(H1/p.dH,p.n_eps))+1;
if F_t~=model.Geometry.NumFaces
    stop
end
factor=1.01;
mesh = generateMesh(model,'Hmax',5*p.dR,'Hmin',p.dmesh*factor,'Hface',{2:F_t,p.dmesh*factor},'GeometricOrder','linear');

% 
% if p.QS==1
%     mesh = generateMesh(model,'Hmax',5*p.dR,'Hface',{2:F_t,p.dmesh*1.01},'GeometricOrder','linear');
% else
%     mesh = generateMesh(model,'Hmax',5*p.dR,'Hface',{2:F_t,p.dmesh*1.01});
% end
end