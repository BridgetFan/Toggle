clear;

% factor=10;
% a1=2.6:0.01:2.7;
% tg=1:1:5;
% a2=1:1:5;
% ret=[0,0.5,0.8,1];

u=[30, 60];
l3=0.1:0.1:0.3;
l4=0.5:0.1:0.8;


[param2,param3]=meshgrid(l3,l4);
param1=ones(size(param2))*u(1);
Date='Jun13_';

parfor i = 1:length(param1(:))
    p1=param1(i);
    p2=param2(i);
    p3=param3(i);
    toggle_fun(Date,i,p1,p2,p3);
end
