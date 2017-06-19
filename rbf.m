%% 2d RBF generator

function [A,B] = rbf(ep,x,y,xi,yi)
%TODO: check xi and yi are the same length

phi = @(ep,r) exp(-ep^2*r.^2);
phi_del = @(ep,r) 4*ep^2*(ep^2*r.^2-1).*exp(-ep^2*r.^2);

[X1,X2] = meshgrid(xi);
[Y1,Y2] = meshgrid(yi);
A = phi(ep,sqrt((X2-X1).^2 + (Y2-Y1).^2));

[X1i,X2i] = meshgrid(xi,x);
[Y1i,Y2i] = meshgrid(yi,y);
B = phi_del(ep,sqrt((X2i-X1i).^2+(Y2i-Y1i).^2));

end

