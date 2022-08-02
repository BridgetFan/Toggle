function c = c_function(location,~,p)
% Specifying the c coefficients as the following vectors.
% For NQS toggle:
%  c1 = [p.D1;p.D1; 0; 0; 0; 0; 0; 0];
%  c2 = [p.D2;p.D2; 0; 0; 0; 0; 0; 0];
% For QS toggle:
%  c1 = [p.D1;p.D1; 0; 0; 0; 0; 0; 0;p.D1_C14;p.D1_C14;p.D1_C4;p.D1_C4];
%  c2 = [p.D2;p.D2; 0; 0; 0; 0; 0; 0;p.D2_C14;p.D2_C14;p.D2_C4;p.D2_C4];

loc_x=round(location.x,5);
loc_y=round(location.y,5);

nr = length(loc_x); % Number of columns
indx=logical((loc_y>-p.Eps).*(loc_x<p.R+p.U*p.dR+p.Eps));  % Locate colony nodes


if p.N==1
    % 1D PDE (NQS w/ Splitting)
    c=zeros(1,nr);
    c(1,:) = p.D2;
    c(1,indx) = p.D1;
else
    c = zeros(2*p.N,nr); % Allocate c
    
    % aTc Diffusion Coeff.
    c(1:2,:) = p.D2;
    c(1:2,indx) = p.D1;
    if p.QS==1
        % QS Signals Diffusion Coeff.
        c((2*p.N-3):(2*p.N-2),:) = p.D2_C14;
        c((2*p.N-3):(2*p.N-2),indx) = p.D1_C14;
        c((2*p.N-1):2*p.N,:) = p.D2_C4;
        c((2*p.N-1):2*p.N,indx) = p.D1_C4;
    end
end

end