function f = f_function(location,state,p)
% Reaction Function
loc_x=round(location.x,5);
loc_y=round(location.y,5);

%N Number of equations
nr = length(loc_x); % Number of columns
f = zeros(p.N,nr); % Allocate f

f(1,:) = p.D2./loc_x.*state.ux(1,:);

if p.QS
    f(p.N-1,:) = p.D2_C14./loc_x.*state.ux(p.N-1,:);
    f(p.N,:) = p.D2_C4./loc_x.*state.ux(p.N,:);
end

indx=logical((loc_y>-p.Eps).*(loc_x<p.R+p.U*p.dR+p.Eps));

dilution=zeros(1,nr);

if sum(indx)
    %plot(loc_x,loc_y,'*');
    u = state.u(1,:);
    x = state.u(2,:);
    y = state.u(3,:);
    yt= state.u(4,:);
    ybar = yt-y;

%     if min([u,x(indx),y(indx),yt(indx),ybar(indx)])<-0.001
%         pause
%     end

    % Now the particular functional form of f
    sqst=-p.kp*u(indx).*y(indx)+p.km*ybar(indx);
    f(1,indx) = p.D1./loc_x(indx).*state.ux(1,indx)+sqst+p.ret*(p.g2+dilution(1,indx)).*ybar(indx);

    if ~p.split
        if p.QS==0
            prod_x=p.A1./(1+(y(indx)/p.theta_y).^2);
            prod_y=p.A2./(1+(x(indx)/p.theta_x).^4);
        elseif p.QS==1
            g=state.u(5,indx);
            h=state.u(6,indx);
            prod_x=p.A1*p.leak1+p.A1.*(1-p.leak1)./(1+(y(indx)/p.theta_y).^2).*(g./(p.theta_g+g));
            prod_y=p.A2./(1+(x(indx)/p.theta_x).^4).*h./(p.theta_h+h);
        end
        f(2,indx) = prod_x;
        f(3,indx) = prod_y+sqst;
        f(4,indx) = prod_y;
    end

    if p.QS
        f(p.N-1,indx) = p.D1_C14./loc_x(indx).*state.ux(p.N-1,indx)+p.A3.*(p.leak3+(1-p.leak3)./(1+(y(indx)/p.theta_y).^2));
        f(p.N,indx)   = p.D1_C4./loc_x(indx).*state.ux(p.N,indx)   +p.A4.*(p.leak4+(1-p.leak4)./(1+(x(indx)/p.theta_x).^4));
    end
end

indx_agar=loc_y<-10*p.dmesh;

% if ~p.split  && sum(indx_agar)
%     if max(max(state.u(2:4,indx_agar)))>2
%         pause
%     end
% end

end

% AGL=logical((loc_y<p.H_agl+p.Eps).*indx);
% if max(AGL)
%     dilution(:,AGL)=p.beta2;
%     AGL_fast=logical(AGL.*(loc_x>p.R+(p.U-1)*p.dR-p.Eps));
%     % plot(loc_x(AGL),loc_y(AGL),'go');
%     if max(AGL_fast)
%         dilution(:,AGL_fast)=p.beta2*log(4)/log(2);
%         %plot(loc_x(AGL_fast),loc_y(AGL_fast),'y*');
%     end
% end


%     if p.split==0
%         x =state.u(2,:);
%         y =state.u(3,:);
%         yt=state.u(4,:);
%     else
%         x = zeros(size(u));y=x;yt=x;
% %         tic
%         x(indx) =griddata(p.mesh.Nodes(1,p.indx),p.mesh.Nodes(2,p.indx),p.ic(p.indx,2),loc_x(indx),loc_y(indx),'nearest');
%         y(indx)=griddata(p.mesh.Nodes(1,p.indx),p.mesh.Nodes(2,p.indx),p.ic(p.indx,3),loc_x(indx),loc_y(indx),'nearest');
%         yt(indx)=griddata(p.mesh.Nodes(1,p.indx),p.mesh.Nodes(2,p.indx),p.ic(p.indx,4),loc_x(indx),loc_y(indx),'nearest');
% %         toc
% 
% %         tic
% %         interpolateSolution(p.ModelSoln_xy,loc_x,loc_y,1:3);
% %         toc
%         % Interp entire thing cost 0.015s; Interp only colony cost 0.005s
%     end
