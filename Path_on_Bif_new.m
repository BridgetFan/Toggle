clear;
filePattern = fullfile('./ClusterResult/QS_FinalVersion/', 'Jun28_3*Soln.mat');
% filePattern = fullfile('Jun21*Soln.mat');
% filePattern = fullfile('./Result/', 'Jul05_17*Soln.mat');
files = dir(filePattern);
for k=1:length(files)
    close all
    clearvars -except filePattern files k
    FileName=files(k).name;
    FileFolder=files(k).folder;
    load([FileFolder '/' FileName]);
    N=p.NT;
    p.N_total=23;
    N0=p.gI*p.T/p.dt;

    ini=nan(p.N,p.T/p.dt*p.N_total+N0+1);

    t=0:length(ini)-1;

    a1=p.A1;
    a2=p.A2;

    p.g1=repmat(p.g1(1),1,length(t));
    p.g2=p.g1;
    % g1=[p.g1;p.g1+p.beta2;p.g1+2*p.beta2];
    p.finer=1;
    theta_x=p.theta_x;
    theta_y=p.theta_y;
    kp=p.kp;
    km=p.km;
    % figure(1),
    % hs=Bif_plot(soln1,p);
    % delete(hs)
    % hold on
    %%

    for ti=1:N
        if exist(['mesh' num2str(ti)],"var")
            for j=1:p.T/p.dt
                %             eval(['p.mesh=mesh' num2str(ti) ';']);
                eval(['p.mesh=mesh' num2str(ti) '{' num2str(j) '};']);
                eval(['soln=soln' num2str(ti) '{' num2str(j) '};']);
                ind=j+(ti-1)*p.T/p.dt+N0+1;

                %N_AGL=findNodes(p.mesh,'region','Face',2);
                H_colony=max(p.mesh.Nodes(2,:));
                R_colony=p.k*H_colony;%max(p.mesh.Nodes(1,N_AGL));

                loc1=find(p.mesh.Nodes(2,:)==H_colony);
                loc2=find((abs(p.mesh.Nodes(2,:))< p.Eps)&(abs(p.mesh.Nodes(1,:)-R_colony)<p.Eps));
                loc=[loc1,loc2];
                soln_u=soln(:,1);
                ini(1:2,ind)=soln_u(loc);
                if p.QS
                    soln_C14=soln(:,end-1);
                    soln_C4=soln(:,end);
                    ini(3:4,ind)=soln_C14(loc);
                    ini(5:6,ind)=soln_C4(loc);
                end
                %         if isfield(p, 'pulse_noisy') && p.AddNoise==1
                %             p.g1(ind)=p.pulse_noisy((ti-1+p.gI)*p.T/p.dt+(0:p.T/p.dt));
                %         end
            end
                            size(ini)

%         else
%             ti=ti+1;
%             break
        end
    end
    p.g2=p.g1;
    u_all_top=ini(1,:);
    u_all_tip=ini(2,:);


    if p.QS
        C14_top=ini(3,:);
        C14_tip=ini(4,:);
        C4_top=ini(5,:);
        C4_tip=ini(6,:);
        C14_top_eff=C14_top./(p.theta_g+C14_top);
        C14_tip_eff=C14_tip./(p.theta_g+C14_tip);
        C4_top_eff=C4_top./(p.theta_h+C4_top);
        C4_tip_eff=C4_tip./(p.theta_h+C4_tip);
    else
        C14_top_eff=1;%ones(size(u_all_top));
        C14_tip_eff=1;%ones(size(u_all_top));
        C4_top_eff=1;%ones(size(u_all_top));
        C4_tip_eff=1;%ones(size(u_all_top));
    end


    g1=[p.g1;p.g1+2*p.beta2];
    g2=[p.g2;p.g2+2*p.beta2];

    A1_top=a1./g1(1,:)./theta_x.*C14_top_eff;
    A2_top=a2./g2(1,:)./theta_y./(kp*u_all_top./(km+g2(1,:))+1).*C4_top_eff;

    if p.gI==1
        ind1=1:p.T*p.dT+1;
        A1_top(ind1)=A1_top(ind1).*g1(1,ind1)./g1(2,ind1);
        A2_top(ind1)=A2_top(ind1).*g1(1,ind1)./g1(2,ind1);
    end

    A1_tip=a1./g1(2,:)/theta_x.*C14_tip_eff;
    A2_tip=a2./g2(2,:)/theta_y./(kp*u_all_tip./(km+g2(2,:))+1).*C4_tip_eff;
    %%

    step=1;
    ind2=[1 step+1:step:length(A1_top)];
    l=length(A1_top(ind2));



    figure,
    plot(p.bif_a1,p.bif_a2,'k-','LineWidth',2);
    hold on
    scatter(A1_top(ind2),A2_top(ind2),100,1:l,'o','filled','DisplayName','Top');
    scatter(A1_tip(ind2),A2_tip(ind2),100,1:l,'^','filled','DisplayName','Rim');
    legend('Location','northwest');
    % colormap parula
    % colormap(flipud(autumn(l-1)))
    % colormap(flipud(copper(l-1)))
    % colormap(flipud(hot(l-1)))
    % colormap(flipud(parula(l-1)))
    colormap(flipud(jet(l-1)))
    caxis([1 l]);
    cbh = colorbar ;
    cbh.Ticks = [1 (l+1)/2 l];
    cbh.TickLabels = sprintfc('%d hr', 0:30:60);
    axis equal;

    a_max=max([A1_top(ind2(end)),A2_top(ind2(end))]);

    % if a_max<3
    %     text(2.3,0.5,'Yellow','HorizontalAlignment','center','FontSize',30);
    %     text(0.5,2.3,'Blue','HorizontalAlignment','center','FontSize',30);
    %     text(2.3,2.3,'Bistable','HorizontalAlignment','center','FontSize',30);
    %     set(gca,'FontSize',30,'XTick',0:3:6,'YTick',0:3:6);
    %     xlim([0 3]);ylim([0 3]);
    % elseif a_max<5
    %     text(3,1,'Yellow','HorizontalAlignment','center','FontSize',30);
    %     text(1,3,'Blue','HorizontalAlignment','center','FontSize',30);
    %     text(3,3,'Bistable','HorizontalAlignment','center','FontSize',30);
    %     title('NQS Toggle');
    %     set(gca,'FontSize',30,'XTick',[0,5],'YTick',[0,5]);
    %     xlim([0 5]);ylim([0 5]);
    % else
    text(4,1,'Yellow','HorizontalAlignment','center','FontSize',30);
    text(1,4,'Blue','HorizontalAlignment','center','FontSize',30);
    text(4,4,'Bistable','HorizontalAlignment','center','FontSize',30);
    set(gca,'FontSize',30,'XTick',0:3:6,'YTick',0:3:6);
    xlim([0 6]);ylim([0 6]);
    % end

    if p.QS
        title('QS Toggle');
    else
        title('NQS Toggle');
    end


    set(gca,'FontSize',30,'XTick',0:5:10,'YTick',0:5:10);
    sc=get(0, 'MonitorPositions');
    po=sc(end,:);
    %     set(gcf,'Position',[po(1) po(2)+100 800 800], 'color', 'w');
    set(gcf,'Position',[po(1) po(2) 800 800], 'color', 'w');
    hold off
    saveas(gcf,[FileFolder '/' p.RunName '_path.png']);
end
% figure,
% plot(u_all_tip,'DisplayName','PDE Simulation Tip');
% hold on
% plot(u_all_top,'DisplayName','PDE Simulation Top');
% % t=1:length(u_all_tip);
% % plot(t,30*exp(-log(2)/20/60*(t-1)),'DisplayName','aTc Degradation');
% plot(t,p.u0*exp(-p.beta1/10*(t)),'DisplayName','aTc Degradation');
%
% legend('Location','northeast');
% set(gcf,'Position',[po(1)+100 po(2)+100 800 800], 'color', 'w');
% set(gca,'FontSize',30,'XTick',0:150*4:150*32,'XTickLabel',num2str((0:10:80)'));
% set(gca,'LineWidth',2)
%
% figure,
% plot(t(1:ind(end))/60,C14_top_eff(1:ind(end)),'-*','Color',p.yellow/2);
% hold on;
% plot(t(1:ind(end))/60,C4_top_eff(1:ind(end)),'-*','Color',p.blue);
% legend('C14 eff','C4 eff','Location','northwest');
% xlim([0, 25]);ylim([0, 1]);
% % xlim([0, t(ind(end))/60]);
% xlabel('Time (hr)');
% set(gca,'FontSize',18);
% stop
%%
clearvars B
sc=get(0, 'MonitorPositions');
p.po=sc(end,:);
% ti_range=1:23;
ti_range=[12,18,24]-p.gI;
st_time=['30';'45';'60'];
% ti_range=17;
rec_dt=1;
frame_i=1;
%close all;
for jj=3:length(ti_range)

    figure(2),
    ti=ti_range(jj);
    % ti=23;
    if exist(['mesh' num2str(ti)],"var")
        for tj=150:150:p.T/p.dt
            p.U=ti-1+floor(round((tj-1)/(p.dt_mesh/p.dt),5))*p.dU;
            eval(['p.mesh=mesh' num2str(ti) '{' num2str(tj) '};']);
            eval(['soln=soln' num2str(ti) '{' num2str(tj) '};']);
            model=createpde(p.N);
            geometryFromMesh(model,p.mesh.Nodes,p.mesh.Elements);
            createPDEResults(model,reshape(soln,[],1));

            time_real=(ti+p.gI-1)*p.T*p.dT+tj*p.dT*p.dt;
            sgt = sgtitle(['t = ' num2str(floor(time_real/60),'%02d') ' hr ' num2str(mod(time_real,60),'%02d') ' min'],'FontSize',20);
            set(sgt,'FontSize', 25,'FontWeight','bold');
            colonyplot(model,soln,p);
            set(gcf,'units','pixels','Position',[p.po(1) p.po(2)+200 1600 1000], 'color', 'w','Resize','off');
                print('-vector','-dsvg',['NQS_aTc' num2str(p.u0) '_t' st_time(jj,:) '_rendered.svg']);
%             pause(1);
            B(frame_i)=getframe(gcf);
            frame_i=frame_i+1;
        end

    else
        jj=jj-1;
        break
    end
end
close all
v = VideoWriter([FileFolder '/' FileName(1:8) '_B'],'MPEG-4');
v.FrameRate = 8;
open(v);
writeVideo(v,B);
close(v);