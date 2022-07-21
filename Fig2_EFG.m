clear;
close all;
a1=3;a2=3;k=1;
load('p_data.mat');
sc=get(0, 'MonitorPositions');
po=sc(end,:);

figure,
scatter(a1/(1+k),a2,100,'^','filled');
hold on;
scatter(a1,a2,100,'o','filled');
scatter(a1,a2/(1+k),100,'s','filled');

plot(p.bif_a1,p.bif_a2,'k-','LineWidth',2);

legend('v0=1, u0=0','v0=0, u0=0','v0=0, u0=1','Cusp Curve');
xlim([0 5]);ylim([0 5]);set(gca,'FontSize',30,'XTick',[0,5],'YTick',[0,5]);
set(gcf,'Position',[po(1) po(2)+200 800 800], 'Resize','off');
box on;

yellow=[201,201,0]/255;
blue=[59,59,252]/255;
%%
u_range=[0 k];
v = 0;
for i=1:length(u_range)
    u=u_range(i);
    figure,
    [X, Y] = meshgrid(0:0.5:5, 0:0.5:5);
    S = a1./(1+(Y/(1+u)).^2)-X; % change this line here
    L = a2./(1+(X/(1+u)).^4)-Y;
    quiver(X, Y, S./sqrt(S.^2+L.^2), L./sqrt(S.^2+L.^2),0.4,'Color',[0 0 0]+0.5,'LineWidth',2)
    hold on
    inp=0:0.01:5;
    plot(a1./(1+(inp./(1+u)).^2),inp,'LineWidth', 5,'Color',yellow);
    plot(inp,a2./(1+(inp./(1+v)).^4),'LineWidth', 5,'Color',blue);
    set(gca,'FontSize',30,'XTick',0:1:5,'YTick',0:1:5);
    axis square;xlim([0 5]);ylim([0 5]);
    xlabel('tilde x');
    ylabel('tilde y');
    legend('Direction Field','x-nullcline','y-nullcline');
    axis tight; xlabel('x'), ylabel('y');
    title(['v0=' num2str(v) ', u0=' num2str(u)]);
    set(gcf,'Position',[po(1) po(2)+200 800 800],'Resize','off');
end