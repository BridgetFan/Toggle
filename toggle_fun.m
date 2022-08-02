%%
% Simulation of the NQS/QS toggle colony coupled ODE-PDE model
% For any questions, please contact 
% Gaoyang (Bridget) Fan at bridget.gfan@gmail.com

%%
close all;
warning('off');
% sympref('HeavisideAtOrigin','default');
tStart = tic;

Date='XX01_';
run_num=1;
p.RunName=[Date num2str(run_num)];

p.QS = 0;                  % Simulate NQS or QS Toggle 

%% Parameters - Space & Time
p.H_agl = 0.01;            % Height of AGL
p.dH = p.H_agl;            % Height growth per time T
p.dmesh = p.dH/2;          % Unit mesh size

p.split = 1;               % Split the ODE var. from the PDE ones

p.T = 150;                 % AGL doubling time (min)
p.dt = 1;                  % PDE-ODE splitting update period (min)
p.dt_mesh = 15*p.dt;       % Colony shape and mesh update period (min)

p.dH_mesh = p.dH*p.dt_mesh/p.T; % Height increment per mesh update dt_mesh

p.dU = p.dt_mesh/p.T;      % Cycle # increment per mesh update dt_mesh

p.gI = 1;                  % Initial cycles #
p.N_total = 24;            % End cycles #
p.NT = p.N_total-p.gI;     % Number of cycles to simulate

p.k = 10;                          % Aspect ratio of the colony (R:H) 
p.dR = p.k*p.dH;

p.Ra = 5;                          % Agar radius (mm)
p.Ha = 3;                          % Agar depth (mm)

p.R_max = p.N_total*p.dR;          % Agar radius at the last cycle
p.H_max = p.N_total*p.dH;          % Agar height at the last cycle

p.R = p.dR;                        % Base colony radius (t=0)
p.H = p.dH;                        % Base colony height (t=0)

%% Parameters - Chemicals

p.u0=30;                   % [aTc] at t=0 (nM)
p.ret = 0.8;               % Portion of aTc returning from aTc-TetR complex

if p.QS==0
    p.N =4;                % # of PDE equations
    p.PDEVar = 1:p.N;
    p.A1 = 100;            % NQS max. LacI production rate (nM/min)
    p.A2 = 250;            % NQS max. TetR production rate (nM/min)
    p.kp = 0.06;           % aTc-TetR binding rate (1/min/nM)
    p.km = p.kp/10;        % aTc-TetR unbinding rate (1/min)

    FileName = [p.RunName '_NQS_' num2str(p.A1) '_' num2str(p.A2) '_' num2str(p.ret*100) '_' num2str(p.u0)];
else
    p.N = 6;               % # of PDE equations
    p.PDEVar = 1:p.N;
    p.A1 = 110;            % QS max. LacI production rate (nM/min)
    p.A2 = 210;            % QS max. TetR production rate (nM/min)
    p.kp = 0.006;          % aTc-TetR binding rate (1/min/nM)
    p.km = p.kp/10;        % aTc-TetR unbinding rate (1/min)

    p.A3 = 250;            % QS max. effective C14 production rate (nM/min)
    p.A4 = 8*p.A3;         % QS max. effective  C4 production rate (nM/min)
    p.leak3=0.2;           % C14 base production percentage
    p.leak4=0.5;           %  C4 base production percentage
    p.theta_g = 200;       % C14 activation threshold (nM)
    p.theta_h = 200;       %  C4 activation threshold (nM)

    FileName = [p.RunName '_QS_' num2str(p.A1) '_' num2str(p.A2) '_' num2str(p.ret*100) '_' num2str(round(p.u0)) '_' num2str(p.A3) '_' num2str(p.A4)];
end

p.beta2 = 1/p.T;                   % Dilution rate (1/min)

% LacI & TetR
p.gamma_1 = log(2)/7;              % Tagged degradation time (1/min) 
p.g1 = p.gamma_1;                  % LacI tagged degradation rate
p.g2 = p.gamma_1;                  % TetR tagged degradation rate
p.theta_x = 500;                   % LacI repression threshold (nM)
p.theta_y = 500;                   % TetR repression threshold (nM)

% aTc
p.beta1 = log(2)/(48*60);          % aTc degradation rate (1/min)
p.D1 = 40*60/10^6;                 % aTc diffusion in colony (mm^2/min)
p.D2 = 10*p.D1;                    % aTc diffusion in agar (mm^2/min)

% C14
p.g_C14 = log(2)/(24*60);          % C14 degradation rate (1/min)
p.D1_C14 = 10*60/10^6;             % C14 diffusion in colony (mm^2/min)
p.D2_C14 = 10*p.D1_C14;            % C14 diffusion in agar (mm^2/min)

% C4
p.g_C4 = p.g_C14;                  % C4 degradation rate (1/min)
p.D1_C4 = 100*60/10^6;             % C4 diffusion in colony (mm^2/min)
p.D2_C4 = 10*p.D1_C4;              % C4 diffusion in agar (mm^2/min)

p.t = 0:p.dt:p.N_total*p.T;

p.AddPulse = 0;            
if p.AddPulse
    t_N1 = (8.3+p.gI-1)*p.T;
    t_N2 = (9.23+p.gI-1)*p.T;
    p.pulse = rectangularPulse(t_N1,t_N2,p.t)*(log(2)/5.178-p.gamma_1)+p.gamma_1;
end

p.Eps = 1e-8;
p.n_eps = -log10(p.Eps)-2;
p.er = [-p.Eps,p.Eps];

%% Initialization
x0 = 1000;                         % LacI initial condition (~ max [LacI])
y0 = 0;                            % TetR initial condition
g0 = 0;                            % C14 initial condition
h0 = 0;                            % C4 initial condition

p.U = p.gI;                        % Record initial cycle #

[model, mesh] = MeshGenerator(p);  % Generate the mesh for colony at t=t_0

p.R = p.gI*p.dR + p.R;             % Radius at t_0
p.H = p.gI*p.dH + p.H;             % Height at t_0
p.U = 0;                           % Reset simulation cycle #
p.mesh = mesh;                    

u0_i = p.u0*exp(-p.beta1*p.gI*p.T);% [aTc] at t_0
ic1 = [u0_i, x0, y0, y0];          % Colony
ic2 = [u0_i, 0, 0, 0];             % Agar

if p.QS == 1
    ic1 = [ic1,g0,h0];
    ic2 = [ic2,g0,h0];
end

p.M = length(mesh.Nodes(1,:));
ind_colony = logical((mesh.Nodes(2,:)>-p.Eps).*(mesh.Nodes(1,:)<(p.R+p.U*p.dR+p.Eps)));

ic = repmat(ic2,p.M,1);                             % Size #M-by-N
ic(ind_colony,:) = repmat(ic1,sum(ind_colony),1);   % Set I.C. for colony

%Define Yellow and Blue Colormap
p.yellow = [201,201,0]/255;
p.blue = [59,59,252]/255;

p.yellow_bright = [246,246,0]/255;
p.blue_bright = [14,14,230]/255;

sc = get(0, 'MonitorPositions');
p.po = sc(end,:);
p.po(2) = p.po(2);

if ~isfield(p,'bif_a1')
    [p.bif_a1, p.bif_a2] = Bif(p);
end

p.N2 = p.T/p.dt;                   % # of splitting update per T
soln_j = {};                       % Recording of solution
%% Solving ODE-PDE
tic
for ti = 1:(p.NT*p.N2)
    tj = mod(ti-1,p.N2)+1;
    p.mesh = mesh;
    p.loc_x = round(p.mesh.Nodes(1,:),p.n_eps);
    p.loc_y = round(p.mesh.Nodes(2,:),p.n_eps);
    p.indx = logical((p.loc_y>-p.Eps).*(p.loc_x<p.R+p.U*p.dR+p.Eps));

    % Calculating Dilution Value for Each Node of the Current Mesh
    loc.x = mesh.Nodes(1,:);
    loc.y = mesh.Nodes(2,:);
    p.dilution = a_function(loc,0,p);
    p.dilution = p.dilution(1,:)-p.beta1;

    model.SolverOptions.ReportStatistics = 'off';
    if p.AddPulse
        p.g1 = p.pulse(ti+p.gI*p.N2);
        p.g2 = p.g1;
    end

    % Solving ODE-PDE
    p.ic = ic;

    model_ic = createPDEResults(model,reshape(ic,[],1));
    setInitialConditions(model,model_ic);
    c_fun = @(state,location) c_function(state,location,p);
    a_fun = @(state,location) a_function(state,location,p);
    f_fun = @(state,location) f_function(state,location,p);
    specifyCoefficients(model,'m',0,'d',1,'c',c_fun,'a',a_fun,'f',f_fun);
    tlist = 0:p.dt;
    Y0(2).soln = [];
    
    % Solving PDE/ODE in parallel
    parfor paral_i = 1:2

        if paral_i == 1
            % Solve PDE
            Y0(paral_i).soln = solvepde(model,tlist);
        else
            % Solve ODE for nodes in Colony (vectorization)
            slices=1:ceil(sum(p.indx)/5):sum(p.indx); % Slicing the vector
            if slices(end)~=sum(p.indx)
                slices(6) = sum(p.indx);
            end
            ind_num = find(p.indx);
            for ode_i = 1:5
                indx_i = ind_num(slices(ode_i):(slices(ode_i+1)-(ode_i<5)));
                [t,Y_i] = ode15s(@(t,X) SplitODE(t,X,p,indx_i),[0 p.dt],ic(indx_i,2:4)');
                % ODE Solution
                Y_i = reshape(Y_i,[],length(indx_i));
                mt = length(t);LacI1_i=Y_i(mt,:);TetR1_i=Y_i(2*mt,:);TetR_T1_i=Y_i(3*mt,:);
                ode_soln_i = [LacI1_i;TetR1_i;TetR_T1_i]';
                Y0(paral_i).soln = [Y0(paral_i).soln;ode_soln_i];
            end
        end
    end

    results = Y0(1).soln;
    X0 = Y0(2).soln;

    % Record PDE Solution
    if p.N>1
        soln_u = results.NodalSolution(:,:,end);
    elseif p.N == 1
        soln_u = results.NodalSolution(:,end);
    end

    % Record both PDE & ODE solutions
    soln = soln_u;
    soln(:,2:4) = 0;
    soln(p.indx,2:4) = X0;
    soln_j{tj} = soln;
    mesh_j{tj} = mesh;

    % I.C. for next update
    ic = soln;

    % Write solution every T
    if ~mod(ti,p.N2)
        eval(['soln' num2str(ti/p.N2) '=soln_j;']);
        eval(['mesh' num2str(ti/p.N2) '=mesh_j;']);%%%%  New
        run_ti = toc;
        disp(['Run' num2str(run_num) ' Cylce' num2str(ti/p.N2) ' took ' num2str(round(run_ti/60,1)) ' min']);
        %%
        if ti == p.N2
            savetofile(tlist,['./Result/' FileName '_Soln.mat'],'tlist',0);
        end
        savetofile(p,['./Result/' FileName '_Soln.mat'],'p',1);  %%%%  New
        savetofile(eval(['soln' num2str(ti/p.N2)]),['./Result/' FileName '_Soln.mat'],['soln' num2str(ti/p.N2)],1);
        savetofile(eval(['mesh' num2str(ti/p.N2)]),['./Result/' FileName '_Soln.mat'],['mesh' num2str(ti/p.N2)],1);
        vars = {['mesh' num2str(ti/p.N2)],['soln' num2str(ti/p.N2)]};%%%%  New
        clear(vars{:});
        soln_j = {};


        figure(1);
        time_real=ti*p.dt+p.gI*p.T;
        set(gcf,'Position',[p.po(1) p.po(2)+200 min(1600,p.po(3)) min(800,p.po(4))], 'color', 'w','Resize','off');
        sgtitle(['t = ' num2str(floor(time_real/60)) ' hr ' num2str(mod(time_real,60)) ' min'],'FontSize',20);
        colonyplot(model,soln,p);
        pause(1);

        tic
    end

    % Update colony shape/mesh at every dt_mesh
    if ~mod(ti,p.dt_mesh/p.dt)
        mesh_old=mesh;
        p.U=round(p.U+p.dU,p.n_eps);
        [model, mesh]=MeshGenerator(p);
        ic=Initialization(mesh_old,mesh,soln,p);
    end
end
tEnd = toc(tStart);

