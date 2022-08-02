function a = a_function(location,~,p)
% Specifying the a coefficients, in our case, the effective degradation
% rate

nr = length(location.x);   % Number of columns
a = zeros(p.N,nr);         % Allocate a

a(1,:)=p.beta1;            % Set aTc degradation rate

if p.QS==1                 % Set AHL degradation rates
    a(p.N-1,:)=p.g_C14;
    a(p.N,:)=p.g_C4;
end

indx=logical((location.y>-p.Eps).*(location.x<p.R+p.U*p.dR+p.Eps));

if ~p.split               % Setting the LacI/TetR degrdation rates 
    a(2,indx)=p.g1;
    a(3,indx)=p.g2;
    a(4,indx)=p.g2;
end

dilution=zeros(p.N,nr);   % Initialing dilution matrix
AGL=logical((location.y<p.H_agl+p.Eps).*indx);  % Locate nodes at AGL layer
if max(AGL)
    dilution(:,AGL)=p.beta2;                    % Set dilution rate at AGL
    AGL_fast=logical(AGL.*(location.x>p.R+(p.U-1)*p.dR-p.Eps)); % Locate rim
    if max(AGL_fast)
        dilution(:,AGL_fast)=2*p.beta2;         % Set dilution rate at rim
    end
end

if p.split
    dilution(2:4,:)=0;    % ODE var. set as constant if splitting
end
a=a+dilution;             % Set a as degradation+dilution
end