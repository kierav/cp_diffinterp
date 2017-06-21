function [u,v] = IC(xx,yy)

%pert
%ss = xx + 0.5*yy;
%tt = xx - 0.5*yy;
%z = (cos(5*th)+2) .* exp(1i*th);
[th2,r2] = cart2pol(xx,yy);
th3 = th2 + pi/8;
z2 = ((cos(3*th3)+5).*exp(-(r2).^2));
%z2 = ((cos(3*th3)+5).* (r2 <= .8));   % .8 good for demo!!
%z2 = ((cos(3*th3)+2) <= 0.5) .* exp(1i*th3);
pert = z2 / max(max(max(z2)));
%max(max(pert))
%min(min(pert))
%pert = exp( -(2*(yy-0.5).^2 + (xx-0.25).^2)  );
%pert = exp(-( 2*(ss-0.5).^2 + 2*(tt-0.7).^2) ) + ...
%       exp(-( 2*(ss).^2 + 2*(tt).^2)  ) + ...
%       exp(-( 2*(ss+0.4).^2 + 2*(tt-.9).^2)  );
%pert = exp( - ( max(abs(ss),abs(tt)) ) );
%pert = pert / max(max(max(pert)));
u0 = 1 - pert;
v0 = 0 + 0.5*pert;

u = u0(:);
v = v0(:);

%[N,temp] = size(u0);
%i = 0;
%plot2d
