%% 2d RBF generator

function [A,B, Ae] = rbf3d(ep,x,y,z,xi,yi,zi)
%x,y : collocation points ; xi,yi : centers
global RBFtype; 
global RBFscale
global RBFpar; 
RBFtype='mq' ;%see frbf.m
RBFpar=-1;%7-1/2;%it works for 'ms'  % Matern/Sobolev (msorder-2/2;)
[X1,X2] = meshgrid(xi);
[Y1,Y2] = meshgrid(yi);
[Z1,Z2] = meshgrid(zi);

[X1i,X2i] = meshgrid(xi,x);
[Y1i,Y2i] = meshgrid(yi,y);
[Z1i,Z2i] = meshgrid(zi,z);

Xcol=[x(:),y(:),z(:)];
Xctr=[xi(:),yi(:),zi(:)];
RBFscale= rand(length(Xctr),1)*0+ep;%
%%Guassian by K
% phi = @(ep,r) exp(-ep.^2*r.^2);
% phi_del = @(ep,r) 4.*ep.^2.*(ep.^2*r.^2-1).*exp(-ep.^2*r.^2);
% 
% 
% A = phi(ep,sqrt((X2-X1).^2 + (Y2-Y1).^2));
% 
% 
% B = phi_del(ep,sqrt((X2i-X1i).^2+(Y2i-Y1i).^2));

%% Meng's code
A=matrixgen( Xctr, Xctr,[0,0,0,0,1], RBFscale);
B=matrixgen( Xcol, Xctr,[1,0,0,0,0], RBFscale);
Ae = matrixgen( Xcol, Xctr, [0,0,0,0,1], RBFscale);

end

