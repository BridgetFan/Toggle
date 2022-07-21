function a = a_function(location,~,p)
% Degradation Coeff.
nr = length(location.x); % Number of columns
a = zeros(p.N,nr); % Allocate a

a(1,:)=p.beta1;

if p.QS==1
    a(p.N-1,:)=p.g_C14;
    a(p.N,:)=p.g_C4;
end

indx=logical((location.y>-p.Eps).*(location.x<p.R+p.U*p.dR+p.Eps));

if ~p.split
    a(2,indx)=p.g1;
    a(3,indx)=p.g2;
    a(4,indx)=p.g2;
end

% figure(2),
% plot(location.x,location.y,'*');

dilution=zeros(p.N,nr);
AGL=logical((location.y<p.H_agl+p.Eps).*indx);
if max(AGL)
    dilution(:,AGL)=p.beta2;
    AGL_fast=logical(AGL.*(location.x>p.R+(p.U-1)*p.dR-p.Eps));
    % plot(location.x(AGL),location.y(AGL),'go');
    if max(AGL_fast)
        dilution(:,AGL_fast)=p.beta2*log(4)/log(2);
        %plot(location.x(AGL_fast),location.y(AGL_fast),'y*');
    end
end
if p.split
    dilution(2:4,:)=0;
end
a=a+dilution;
end