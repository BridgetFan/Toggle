function f = f_function(location,state,p)
% Specifying the f coefficients, the reaction functions

% Read the location of the Gauss points
loc_x=round(location.x,5);
loc_y=round(location.y,5);

% N Number of equations
nr = length(loc_x); % Number of columns
f = zeros(p.N,nr);  % Allocate f

% Adjustment term to account for radial symmetry reduction 
f(1,:) = p.D2./loc_x.*state.ux(1,:);                     % aTc
if p.QS
    f(p.N-1,:) = p.D2_C14./loc_x.*state.ux(p.N-1,:);     % C14
    f(p.N,:) = p.D2_C4./loc_x.*state.ux(p.N,:);          % C4
end

% Locate points in the colony
indx=logical((loc_y>-p.Eps).*(loc_x<p.R+p.U*p.dR+p.Eps));

% Reactions in the colony subdomain
if sum(indx)
    u = state.u(1,:);
    x = state.u(2,:);
    y = state.u(3,:);
    yt= state.u(4,:);
    ybar = yt-y;

    
    sqst=-p.kp*u(indx).*y(indx)+p.km*ybar(indx);     % aTc-TetR (un)binding
    f(1,indx)=p.D1./loc_x(indx).*state.ux(1,indx)+sqst+p.ret*p.g2.*ybar(indx);

    % If no splitting PDE & ODE
    if ~p.split
        if p.QS==0
            prod_x=p.A1./(1+(y(indx)/p.theta_y).^2);
            prod_y=p.A2./(1+(x(indx)/p.theta_x).^4);
        elseif p.QS==1
            g=state.u(5,indx);
            h=state.u(6,indx);
            prod_x=p.A1./(1+(y(indx)/p.theta_y).^2).*(g./(p.theta_g+g));
            prod_y=p.A2./(1+(x(indx)/p.theta_x).^4).*(h./(p.theta_h+h));
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

end
