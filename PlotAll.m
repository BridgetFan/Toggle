clear;
close all;
tic
filePattern = fullfile('./ClusterResult', 'Jun09_1*Soln.mat');
% filePattern = fullfile('./Result', 'May25_25*Soln.mat');
% filePattern = fullfile('./Result', 'Mar20*30*.mat');

files = dir(filePattern);
rec_dt=3;
frame_i=1;
sc=get(0, 'MonitorPositions');
p.po=sc(1,:);

%% Plot and Recording
for k=1:length(files)
    FileName=files(k).name;
    FileFolder=files(k).folder;
    load([FileFolder '/' FileName]);
    sc=get(0, 'MonitorPositions');
    p.po=sc(end,:);

    if ~isfield(p,'leak2')
        p.leak2=0;
    end
    if ~isfield(p,'bif_a1')
        [p.bif_a1, p.bif_a2]=Bif(p);
        save([FileFolder '/' FileName],'p','-append'); 
    end

    % Plotting from Cycle ti_0 to ti_1
    ti_0=1;
    ti_1=p.NT;

    % In each cycle, plot every rec_dt unit time
    i_0=1;
    i_1=p.T/p.finer/rec_dt+1;

    for ti=ti_0:ti_1
        p.U=ti-1;
        eval(['p.mesh=mesh' num2str(ti) ';']);
        eval(['soln=soln' num2str(ti) ';']);
        model=createpde(p.N);
        geometryFromMesh(model,p.mesh.Nodes,p.mesh.Elements);
        createPDEResults(model,reshape(soln(:,1:p.N,1),[],1));
        for i=i_0:i_1
            ind=round(1+(i-1)*rec_dt/p.dt,5);
            soln_i=soln(:,:,ind);
            time_real=(i-1)*p.dT*rec_dt+(ti-1+p.gI)*p.T/p.finer*p.dT;
            if isfield(p, 'pulse_noisy')
                p.g1 = p.pulse_noisy((ti-1+p.gI)*p.T/p.finer/p.dt+i-1);
                p.g2 = p.g1;
            end
            % time_real=results.SolutionTimes(ind)*p.dT+(p.U+gI)*T*p.dT;
%             figure(1);
%             sgtitle(['t = ' num2str(floor(time_real/60)) ' hr ' num2str(mod(time_real,60)) ' min'],'FontSize',20);
%             if p.QS==0
%                 figplot(model,soln_i,p);
%             elseif p.QS==1
%                 figplot_5(model,soln_i,p);
%             end
%             A(frame_i)=getframe(gcf);

            figure(2);
            set(gcf,'Position',[p.po(1) p.po(2)+200 min(1600,p.po(3)) min(800,p.po(4))], 'color', 'w','Resize','off');
            sgtitle(['t = ' num2str(floor(time_real/60)) ' hr ' num2str(mod(time_real,60)) ' min'],'FontSize',20);
            colonyplot(model,soln_i,p);
            B(frame_i)=getframe(gcf);
            pause(1);
%             pause 

            %             figure(3);
            %             ind1r = findNodes(p.mesh,'box',[p.R+(p.U-1)*p.dR p.R+p.U*p.dR]+p.er,[0 p.dH]+p.er);
            %             dil=p.beta2*log(4)/log(2);
            % %             ind1r = (p.mesh.Nodes(2,:)==max(p.mesh.Nodes(2,:)));dil=0;
            %             plot(soln_i(ind1r,2),soln_i(ind1r,3),'r*');
            %             hold on;
            %             ind_var=0:0.1:2500;
            %
            %             g_eff=p.g1+dil;
            %             x=p.A1*(p.leak1+(1-p.leak1)./(1+(ind_var/p.theta_y).^2))./g_eff;
            %             y_tot=p.A2./(1+(ind_var/p.theta_x).^4)./g_eff;
            %             y_free=y_tot./(p.kp*soln_i(ind1r,1)./(p.km+g_eff)+1);
            %             plot(x,ind_var,'Color',p.yellow);
            %             plot(ind_var,y_free,'Color',p.blue);
            %             hold off
            %             set(gcf,'Position',[1200 400 600 600],'color', 'w');

            frame_i=frame_i+1;
        end
        p.U=p.U+1;
%         save([FileFolder '/' FileName(1:end-9) '_B_Bif'],'B');
        %         if ~mod(ti+p.gI,4)
        %             savefig(figure(1),[FileName(1:end-5) '_A_' num2str((ti+p.gI)*T/6) '.fig']);
        %             savefig(figure(2),[FileName(1:end-5) '_B_' num2str((ti+p.gI)*T/6) '.fig']);
        %         end
    end
    toc
%             v = VideoWriter([FileFolder '/' FileName '_A'],'MPEG-4');
%             v.FrameRate = 5;
%             open(v);
%             writeVideo(v,A);
%             close(v);

    v = VideoWriter([FileFolder '/' FileName(1:end-9) '_B_Bif'],'MPEG-4');
    v.FrameRate = 5;
    open(v);
    writeVideo(v,B);
    close(v);
%     save(['./Result/' FileName '_A.mat'],'A');
%     save(['./Result/' FileName '_B.mat'],'B');
end

% for i=1:23
%     eval(['soln=soln' num2str(i) ';']);
%     figure(1),
%     Bif_plot(soln,p);
%     set(gcf,'Position',[p.po(1) p.po(2) 1000 800], 'color', 'w');
%     pause
% end

%%
%         ind_top=(mesh.Nodes(2,:)==max(mesh.Nodes(2,:)));
%         p.track_top=[p.track_top;soln_i(ind_top,:)];
%         p.track_bot=[p.track_bot;soln_i(size(soln,1),:)];
%         p.track_agar=[p.track_agar;soln_i(ind_agar,:)];

%             saveas(gcf,['Result/' FileName '_' num2str(ti)]);

