% function toggle_fun(Date,run_num,param1,param2,param3)
% clearvars -except p
sympref('HeavisideAtOrigin','default');
close all;
warning('off');
% profile on
tStart = tic;


p.u0=param1;%12.5/425.5*1e3; %param1;
p.leak3=param2;%0.2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.leak4=param3;%0.5;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date='Jun13_';
% run_num=2;
p.finer=1;
p.H_agl=0.01;
p.dH=p.H_agl/p.finer; %unit hight growth
p.dmesh=p.dH/p.finer;
p.split=1;

p.RunName=[Date num2str(run_num)];
p.QS = 1;
p.T=15;
p.dt=0.1;            % splitting update frequency
p.dt_mesh=15*p.dt;   % colony/mesh update frequency

p.dH_mesh=p.dH*p.dt_mesh/p.T;
p.dU=p.dt_mesh/p.T;

p.AddNoise=0;
p.factor=0.6022;

p.dT=10; % unit time (min)
p.gI=1*p.finer;
p.N_total=15*p.finer;%24%%%%%31%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.NT=p.N_total-p.gI;

p.gamma_1=log(2)/7*p.dT;

if p.QS
    p.kp=0.06*p.dT/10;
else
    p.kp=0.06*p.dT/1;
end
p.km=p.kp/10;


p.ret=0.8;%param2;

p.beta2=1/p.T;     %%%%%%% not log(2)

%% Noiseb
p.t=round(0:p.dt:p.N_total*p.T,p.n_eps);

if p.AddNoise
    t_N1=(8.3+p.gI-1)*p.T;
    t_N2=(9.23+p.gI-1)*p.T;
    p.pulse=rectangularPulse(t_N1,t_N2,p.t)*(log(2)/5.178*p.dT-p.gamma_1)+p.gamma_1;
%     p.pulse=rectangularPulse(12*p.T,13*p.T,p.t)*(log(2)/6.2*p.dT-p.gamma_1)+p.pulse;
end


%
% figure,
% plot(p.t*p.dT/60,p.pulse_noisy);
% hold on
% % yline(log(2)/7*p.dT);
% % yline(log(2)/6*p.dT);
% set(gca,'FontSize',20,'XTick',0:10:p.t(end)*p.dT/60,'XTick',0:20:60);
% xlabel('Time (hr)');
% ylim([0.5 1.5]);



if p.QS==1
    p.leak1=0; % Leakiness%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p.leak2=0;
else
    p.leak1=0;
    p.leak2=0;
end


if p.QS==0
    p.N =4;   % # of PDE equations
    p.PDEVar=1:p.N;
    p.A1=100*p.dT;
    p.A2=250*p.dT;
    FileName=[p.RunName '_NQS_' num2str(p.A1/p.dT) '_' num2str(p.A2/p.dT) '_' num2str(p.ret*100) '_' num2str(p.u0)];
else
    %     if p.split==1
    %         p.N=3;
    %         p.PDEVar=[1,5,6];
    %     else
    p.N =6;   % # of PDE equations
    p.PDEVar=1:p.N;
    %     end
    p.A1=110*p.dT/1.2; %%%%% Changed
    p.A2=210*p.dT/1.2; %%%%% Changed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p.A3=250*p.dT; %%%%% Changed
    p.A4=8*p.A3;
    FileName=[p.RunName '_QS_' num2str(p.A1) '_' num2str(p.A2) '_' num2str(p.ret*100) '_' num2str(round(p.u0)) '_' num2str(p.A3) '_' num2str(p.A4)];
    p.theta_g=2e2;
    p.theta_h=2e2;
end

% if AddNoise == 1
%     FileName=[FileName '_' num2str(p.factor)];
% else
%     FileName=[FileName '_0'];
% end




% LacI & TetR
p.g1=p.gamma_1;
p.g2=p.g1;
p.theta_x=500;
p.theta_y=500;

% aTc
p.beta1=log(2)/(48*60)*p.dT;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.D1=40*60*p.dT/10^6;  % colony
p.D2=10*p.D1; % agar

% AHL
p.D1_C14=10*60*p.dT/10^6;% colony %%%%% Changed
p.D2_C14=10*p.D1_C14; % agar %%%%% Changed
p.g_C14=log(2)/(24*60)*p.dT; %%%%% Changed


p.D1_C4=100*60*p.dT/10^6;% colony %%%%% Changed
p.D2_C4=10*p.D1_C4;% agar %%%%% Changed
p.g_C4=p.g_C14;


%% Initial Condition
x0=1000;%300
y0=0;%p.A2/p.theta_y/p.gamma_1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g0_1=0;
h0_1=0;
g0_2=0;
h0_2=0;


p.Eps=1e-8;
p.n_eps=-log10(p.Eps)-2;

%% Domain
p.k=10;  % Aspect Ratio of the Colony


% p.dt=0.2;

p.dR=p.k*p.dH;

p.Ra = 5;  %It has to be p.dH*2^n in order for the adaptive mesh work in the agar
p.Ha = 3;


p.R_max=(p.N_total+p.finer)*p.dR;
p.H_max=(p.N_total+p.finer)*p.dH;


p.R = p.dR;
p.H = p.dH;
p.U = p.gI;

[model, mesh]=MeshGenerator(p);


p.H = p.gI*p.dH + p.H;% + p.U*p.dH;
p.R = p.gI*p.dR + p.R;%
p.U = 0;
p.mesh=mesh;

%% Initial Condition
ic1 = [p.u0*exp(-p.beta1*p.gI*p.T), x0, y0, y0]; % Colony
ic2 = [p.u0*exp(-p.beta1*p.gI*p.T), 0, 0, 0]; % Agar

if p.QS==1
    ic1 = [ic1,g0_1,h0_1];
    ic2 = [ic2,g0_2,h0_2];
end

p.M=length(mesh.Nodes(1,:));
ind_wb = logical((mesh.Nodes(2,:)>-p.Eps).*(mesh.Nodes(1,:)<(p.R+p.U*p.dR+p.Eps)));

ic=repmat(ic2,p.M,1);   % Size #M-by-N
ic(ind_wb,:)=repmat(ic1,sum(ind_wb),1);

%% Plot I.C.
%Define Yellow and Blue Colormap
p.yellow=[201,201,0]/255;
p.blue=[59,59,252]/255;

p.yellow_bright=[246,246,0]/255;
p.blue_bright=[14,14,230]/255;

p.er=[-p.Eps,p.Eps];

sc=get(0, 'MonitorPositions');
p.po=sc(end,:);
p.po(2)=p.po(2);

%% Bif Curve Calculation
if ~isfield(p,'bif_a1')
    [p.bif_a1, p.bif_a2]=Bif(p);
end

p.N2=p.T/p.dt;
% soln_j{1}=ic;
soln_j={};
%%

% t0=floor(t_N1);
% ti_0=(t0-1)*p.N2+1;
% if t0~=1
%     path=pwd;
%     load([path '/Result/Jul05_1_QS_1100_2100_80_45_2500_20000_Soln.mat'], ['mesh' num2str(t0-1)]);
%     load([path '/Result/Jul05_1_QS_1100_2100_80_45_2500_20000_Soln.mat'], ['soln' num2str(t0-1)]);
%     eval(['mesh=mesh' num2str(t0-1) '{150};']);
%     eval(['soln=soln' num2str(t0-1) '{150};']);
%     mesh_old=mesh;
%     p.U=t0-1;
%     [model, mesh]=MeshGenerator(p);
%     ic=Initialization(mesh_old,mesh,soln,p);
%     savetofile(p,['./Result/' FileName '_Soln.mat'],'p',0);
% end
tic
%%
p.NT=23;
ti_0=1;
for ti=ti_0:(p.NT*p.N2)
    tj=mod(ti-1,p.N2)+1;
    p.mesh=mesh;
    p.loc_x=round(p.mesh.Nodes(1,:),p.n_eps);
    p.loc_y=round(p.mesh.Nodes(2,:),p.n_eps);
    p.indx=logical((p.loc_y>-p.Eps).*(p.loc_x<p.R+p.U*p.dR+p.Eps));

    % Calculating Dilution Value for Current Mesh
    loc.x=mesh.Nodes(1,:);
    loc.y=mesh.Nodes(2,:);
    p.dilution=a_function(loc,0,p);
    p.dilution=p.dilution(1,:)-p.beta1;


    model.SolverOptions.ReportStatistics = 'off';
    %     soln_j{1}=ic;
    if p.AddNoise
%         p.g1 = p.pulse(ti+p.gI*p.N2);
        p.g1 = p.pulse(end);
        p.g2 = p.g1;
    end

    % Solve PDE part
    p.ic=ic;

    model_ic = createPDEResults(model,reshape(ic,[],1));
    setInitialConditions(model,model_ic);
    c_fun=@(state,location) c_function(state,location,p);
    a_fun=@(state,location) a_function(state,location,p);
    f_fun=@(state,location) f_function(state,location,p);
    specifyCoefficients(model,'m',0,'d',1,'c',c_fun,'a',a_fun,'f',f_fun);
    tlist = 0:(p.dt/1):p.dt;%%%%%%%%%%%%p.dt/5
    Y0(2).soln = [];

    parfor (paral_i=1:2, 2)
%     for paral_i=1:2
        if paral_i==1
            Y0(paral_i).soln = solvepde(model,tlist);
        else
            % Solve ODE part
            % Run ODE for nodes in Colony
            slices=1:ceil(sum(p.indx)/5):sum(p.indx);
            if slices(end)~=sum(p.indx)
                slices(6)=sum(p.indx);
            end
            ind_num=find(p.indx);
            for ode_i=1:5
                indx_i=ind_num(slices(ode_i):(slices(ode_i+1)-(ode_i<5)));
                [t,Y_i]=ode15s(@(t,X) SplitODE(t,X,p,indx_i),[0 p.dt],ic(indx_i,2:4)');
                % ODE Solution
                Y_i = reshape(Y_i,[],length(indx_i));
                mt = length(t);LacI1_i=Y_i(mt,:);TetR1_i=Y_i(2*mt,:);TetR_T1_i=Y_i(3*mt,:);
                ode_soln_i=[LacI1_i;TetR1_i;TetR_T1_i]';
                Y0(paral_i).soln=[Y0(paral_i).soln;ode_soln_i];
            end
        end
    end

    results=Y0(1).soln;
    X0=Y0(2).soln;%X0=vertcat(Y0(2:6).soln);

    % PDE Solution
    if p.N>1
        soln_u = results.NodalSolution(:,:,end);
    elseif p.N==1
        soln_u = results.NodalSolution(:,end);
    end


    % Initializing
    %             soln=soln_u(:,[1,1,1,1:(p.N-3)]);
    soln=soln_u;
    soln(:,2:4)=0;
    soln(p.indx,2:4)=X0;
    ic=soln;
    soln_j{tj}=soln;
    mesh_j{tj}=mesh;


    if ~mod(ti,p.N2)
        eval(['soln' num2str(ti/p.N2) '=soln_j;']);
        eval(['mesh' num2str(ti/p.N2) '=mesh_j;']);%%%%  New
        run_ti=toc;
        disp(['Run' num2str(run_num) ' Cylce' num2str(ti/p.N2) ' took ' num2str(round(run_ti/60,1)) ' min']);
        %%
        if ti==p.N2
            savetofile(tlist,['./Result/' FileName '_Soln.mat'],'tlist',0);
        end
        savetofile(p,['./Result/' FileName '_Soln.mat'],'p',1);  %%%%  New
        savetofile(eval(['soln' num2str(ti/p.N2)]),['./Result/' FileName '_Soln.mat'],['soln' num2str(ti/p.N2)],1);
        savetofile(eval(['mesh' num2str(ti/p.N2)]),['./Result/' FileName '_Soln.mat'],['mesh' num2str(ti/p.N2)],1);
        vars= {['mesh' num2str(ti/p.N2)],['soln' num2str(ti/p.N2)]};%%%%  New
        clear(vars{:});
        soln_j={};


        time_real=(ti*p.dt+p.gI*p.T)*p.dT;
%         figure(1);
%         sgtitle(['t = ' num2str(floor(time_real/60)) ' hr ' num2str(mod(time_real,60)) ' min'],'FontSize',20);
%         if p.QS==0
%             figplot(model,soln,p);
%         elseif p.QS==1
%             figplot_5(model,soln,p);
%             pause(0.1);
% %             set(gcf,'Position',[p.po(1) p.po(2)+200 min(1200,p.po(3)) min(800,p.po(4))],'color', 'w','Resize','off');
%             set(gcf,'Position',[p.po(1) p.po(2)+50 min(1200,p.po(3)) min(700,p.po(4))],'color', 'w','Resize','off');
%         end


        figure(2);
%         set(gcf,'Position',[p.po(1) p.po(2)+200 min(1600,p.po(3)) min(800,p.po(4))], 'color', 'w','Resize','off');
        set(gcf,'Position',[p.po(1) p.po(2)+50 min(1400,p.po(3)) min(700,p.po(4))], 'color', 'w','Resize','off');
        sgtitle(['t = ' num2str(floor(time_real/60)) ' hr ' num2str(mod(time_real,60)) ' min'],'FontSize',20);
        colonyplot(model,soln,p);
        pause(1);

        tic
    end
    if ~mod(ti,p.dt_mesh/p.dt)
        mesh_old=mesh;
        p.U=round(p.U+p.dU,p.n_eps);
        [model, mesh]=MeshGenerator(p);
        ic=Initialization(mesh_old,mesh,soln,p);
    end


end
tEnd = toc(tStart);
% save([p.RunName '_partdata.mat'],'p','soln_j');
% save(['./Result/' FileName '_Soln.mat'],p,'-append');
% save(['./Result/' p.RunName '_workspace.mat']);

% end


