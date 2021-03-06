clear all; clf
%% parameters and functions of the PDE
F = 0.054;  k = 0.063;  nu = 1/30^2;  eta = nu/3;
f = @(u,v) (-u.*v.*v  +  F*(1-u));
g = @(u,v) ( u.*v.*v  -  (F+k)*v);

useLocal = 1; % 1 - use 4x4 points, 0 - use all points

%% make a grid in x,y and choose a timestep
ep = .75;
N = 80;
dx = 4/N;
x1d = (-2:dx:2-dx)';
y1d = x1d;
dt = 1/2;
[xx, yy] = meshgrid(x1d, y1d);

%% closest point function
[cpx,cpy,dist] = cpBeanCurve(xx, yy);
cpx = cpx(:);
cpy = cpy(:);
% make into vectors
xg = xx(:); yg = yy(:);
cpxg = cpx(:); cpyg = cpy(:);

%% Banding: do calculation in a narrow band around the circle
dim = 2;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band);
xg = xx(band); yg = yy(band);

%% initial condtions
[u, v] = IC(xg,yg);

%% Construct an interpolation matrix for closest point
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');

E = interp2_matrix(x1d, y1d, cpxg, cpyg, p, band);
[Ei,Ej,Es] = interp2_matrix(x1d,y1d,cpxg,cpyg,p, band);
Ej = reshape(Ej,length(cpxg),(p+1)^2);

%L = laplacian_2d_matrix(x1d,y1d, order, band);

% construct rbf matrix
D = sparse(length(cpxg),length(cpxg));
if useLocal == 1
    for j = 1:length(cpxg)
        x = xg(Ej(j,:));
        y = yg(Ej(j,:));
        [A,B] = rbf(ep,cpxg(j),cpyg(j),x,y);
        D(j,Ej(j,:)) = B*pinv(A);
        %D(j,Ej(j,:)) = B/A;
    end
else
    [A,B] = rbf(ep,cpxg,cpyg,xg,yg);
    D = B*pinv(A);   
end

%% time loop
for i=1:10000
  unew = u + dt*f(u,v) + dt*nu*(D*u);
  vnew = v + dt*g(u,v) + dt*eta*(D*v);
  
  % Check any difference on the value after running for a period
  %diff = E*unew - u
  
  u = E*unew;
  v = E*vnew;
  
  t = i*dt;
  
  % plot over computation band
  %if (i < 60) skip = 1; else skip = 20; end
  skip = 25;

  % plot every 'skip' timesteps
  if (i==1) || (mod(i,skip) == 0)
    plot2d_compdomain(u, xg, yg, dx, dx, 1)
    title( ['embedded domain: soln at time ' num2str(t) ...
            ', timestep #' num2str(i)] );
    xlabel('x'); ylabel('y');
    colormap(jet);
    drawnow
    %pause
  end
  %u = Eplot*unew;
  %v = Eplot*vnew;
  
  %make_plots
end
