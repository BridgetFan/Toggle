filePattern = fullfile('./ClusterResult/', 'Jul09_2*Soln.mat');
files = dir(filePattern);
%%
for k=1:length(files)
    FileName=files(k).name;
    FileFolder=files(k).folder;
    load([FileFolder '/' FileName]);
    sc=get(0, 'MonitorPositions');
    p.po=sc(1,:);
    ti_range=1:23;
    rec_dt=1;
    frame_i=1;
    for jj=1:length(ti_range)

        figure(1),
        ti=ti_range(jj);

        if exist(['mesh' num2str(ti)],"var")
            if ti==1
                range=[1 15:15:p.T/p.dt];
            else
                range=15:15:p.T/p.dt;
            end
            for tj=range
                p.U=ti-1+floor(round((tj-1)/(p.dt_mesh/p.dt),5))*p.dU;
                eval(['p.mesh=mesh' num2str(ti) '{' num2str(tj) '};']);

                if tj == 1
                    if ~exist('ic',"var")
                    ind_wb = logical((p.mesh.Nodes(2,:)>-p.Eps).*(p.mesh.Nodes(1,:)<(p.R+p.U*p.dR+p.Eps)));
                    ic1 = [p.u0*exp(-p.beta1*p.gI*p.T), 1000, 0, 0]; % Colony
                    ic2 = [p.u0*exp(-p.beta1*p.gI*p.T), 0, 0, 0]; % Agar
                    if p.QS==1
                        ic1 = [ic1,0,0];
                        ic2 = [ic2,0,0];
                    end
                    ic=repmat(ic2,length(p.mesh.Nodes(1,:)),1);   % Size #M-by-N
                    ic(ind_wb,:)=repmat(ic1,sum(ind_wb),1);
                    end
                    soln=ic;
                    time_real=(ti+p.gI-1)*p.T*p.dT;
                else
                    eval(['soln=soln' num2str(ti) '{' num2str(tj) '};']);
                    time_real=(ti+p.gI-1)*p.T*p.dT+tj*p.dT*p.dt;
                end

                model=createpde(p.N);
                geometryFromMesh(model,p.mesh.Nodes,p.mesh.Elements);
                createPDEResults(model,reshape(soln,[],1));

                sgt = sgtitle(['t = ' num2str(floor(time_real/60),'%02d') ' hr ' num2str(mod(time_real,60),'%02d') ' min'],'FontSize',20);
                set(sgt,'FontSize', 25,'FontWeight','bold');
                colonyplot(model,soln,p);
                set(gcf,'units','pixels','Position',[p.po(1) p.po(2) 1600 1000], 'color', 'w','Resize','off');
                B(frame_i)=getframe(gcf);
                frame_i=frame_i+1;
            end
        end
    end
    save([FileFolder '/' FileName(1:end-4) '_B.mat'],'B','-v7.3');
    close all
    v = VideoWriter([FileFolder '/' FileName(1:end-4) '_B'],'MPEG-4');
    v.FrameRate = 8;
    open(v);
    writeVideo(v,B);
    close(v);
end