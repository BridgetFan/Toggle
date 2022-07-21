clear;
files=dir(fullfile('./Result/','Jun27_7_NQS_*'));
FileName=files.name;
FileFolder=files.folder;
load([FileFolder '/' FileName]);
vars = who(); 
TF1 = contains(vars, 'mesh');
t0=numel(vars(TF1))+1;
ti_0=(t0-1)*p.N2+1;
eval(['mesh=mesh' num2str(t0-1) '{150};']);
eval(['soln=soln' num2str(t0-1) '{150};']);
mesh_old=mesh;
p.U=t0-1;
[model, mesh]=MeshGenerator(p);
ic=Initialization(mesh_old,mesh,soln,p);
% savetofile(p,['./Result/' FileName],'p',0);
TF2 = contains(vars, 'soln'); 
vars=vars(TF1|TF2);
clear(vars{:});

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
        p.g1 = p.pulse(ti+p.gI*p.N2);
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
            savetofile(tlist,['./Result/' FileName],'tlist',0);
        end
        savetofile(p,['./Result/' FileName],'p',1);  %%%%  New
        savetofile(eval(['soln' num2str(ti/p.N2)]),['./Result/' FileName],['soln' num2str(ti/p.N2)],1);
        savetofile(eval(['mesh' num2str(ti/p.N2)]),['./Result/' FileName],['mesh' num2str(ti/p.N2)],1);
        vars= {['mesh' num2str(ti/p.N2)],['soln' num2str(ti/p.N2)]};%%%%  New
        clear(vars{:});
        soln_j={};


%         time_real=(ti*p.dt+p.gI*p.T)*p.dT;
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


%         figure(2);
% %         set(gcf,'Position',[p.po(1) p.po(2)+200 min(1600,p.po(3)) min(800,p.po(4))], 'color', 'w','Resize','off');
%         set(gcf,'Position',[p.po(1) p.po(2)+50 min(1400,p.po(3)) min(700,p.po(4))], 'color', 'w','Resize','off');
%         sgtitle(['t = ' num2str(floor(time_real/60)) ' hr ' num2str(mod(time_real,60)) ' min'],'FontSize',20);
%         colonyplot(model,soln,p);
%         pause(1);

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


