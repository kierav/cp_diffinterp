%% Reaction Diffusion on a sphere
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.
%
% This example solves the heat equation on the surface of a sphere,
% with initial conditions u = cos(4*theta)

cpf = @cpSphere;
paramf = @paramSphere;

useLocal = 1; % 1 - use 4x4x4 points, 0 - use all points

% 3D example on a sphere
% Construct a grid in the embedding space
ep = 1;
dx = 0.1;                   % grid size

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

nx = length(x1d);
ny = length(y1d);
nz = length(z1d);

%% Find closest points on the surface
% For each point (x,y,z), we store the closest point on the sphere
% (cpx,cpy,cpz)

% meshgrid is only needed for finding the closest points, not afterwards
[xx yy zz] = meshgrid(x1d, y1d, z1d);
% function cpSphere for finding the closest points on a sphere
[cpx, cpy, cpz, dist] = cpf(xx,yy,zz);
% make into vectors
cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);


%% Banding: do calculation in a narrow band around the sphere
dim = 3;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
xg = xx(band); yg = yy(band); zg = zz(band);

%% Construct an interpolation matrix for closest point

% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');

E = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, p, band);
[Ei,Ej,Es] = interp3_matrix(x1d,y1d,z1d,cpxg,cpyg,cpzg,p,band);
Ej = reshape(Ej,length(cpxg),(p+1)^3);
% e.g., closest point extension:
%u = E*u;

%% Create Laplacian matrix for heat equation
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);

%% construct RBF matrix
D = sparse(length(cpxg),length(cpxg));
if useLocal == 1
    for j = 1:length(cpxg)
        x = xg(Ej(j,:));
        y = yg(Ej(j,:));
        z = zg(Ej(j,:));
        [A,B] = rbf3d(ep,cpxg(j),cpyg(j),cpzg(j),x,y,z);
        D(j,Ej(j,:)) = B*pinv(A);
%         D(j,Ej(j,:)) = B/A;
    end
else
    [A,B] = rbf3d(ep,cpxg,cpyg,cpzg,xg,yg,zg);
    D = B*pinv(A);   
end

%%plots the stability region of FE
%figure
%plot(real(eig(D)),imag(eig(D)),'*')
%pause 
%hold on
%z = exp(1i*pi*(0:200)/100); r = z-1;
%plot(r/dt, 'r')
%pause
   
%% Construct an interpolation matrix for plotting on sphere

% plotting grid on sphere, based on parameterization
[xp,yp,zp] = paramf(256);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, r] = cart2sph(xp1,yp1,zp1);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);

figure(1); clf;

%% parameters and functions for Gray--Scott
FF = 0.054;  kk = 0.06;  nuu = 1/(3/dx)^2;  nuv = nuu/3;
f = @(u,v) (-u.*v.*v  +  FF*(1-u));
g = @(u,v) ( u.*v.*v  -  (FF+kk)*v);
% nuv = 8.87*10^-4; nuu = 0.516*nuv; 
% alpha = 0.899; beta = -0.91; gamma = -0.899;
% tau1 = 3.5; tau2 = 0;
% f = @(u,v) alpha*u.*(1-tau1*v.*v) + v.*(1-tau2*u);
% g = @(u,v) beta*v.*(1+alpha*tau1/beta*u.*v) + u.*(gamma+tau2*v);

%% initial conditions - small perturbation from steady state
pert = 0.5*exp(-(10*(zg-.1)).^2) + 0.5*rand(size(xg));
u0 = 1 - pert;  v0 = 0.5*pert;
u = u0;  v = v0;
u_rbf = u0; v_rbf = v0;

%% Time-stepping for the heat equation

Tf = 1000;
dt = 1/6*(1/max(nuu,nuv))*dx^2;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;

figure(1);
sphplot = Eplot*u;
sphplot = reshape(sphplot, size(xp));
Hplot = surf(xp, yp, zp, sphplot);
title('initial u')
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
view(-10, 60)
axis off;
shading interp
camlight left
colorbar

for kt = 1:numtimesteps
    % explicit Euler timestepping
    unew = u + dt*(nuu*(L*u) + f(u,v));
    unew_rbf = u_rbf + dt*(nuu*(D*u_rbf) + f(u_rbf,v_rbf));
    vnew = v + dt*(nuv*(L*v) + g(u,v));
    vnew_rbf = v_rbf + dt*(nuv*(D*v_rbf) + g(u_rbf,v_rbf));
    
    % closest point extension
    if ( mod(kt,1) == 0)
        u = E*unew;
        v = E*vnew;
        u_rbf = E*unew_rbf;
        v_rbf = E*vnew_rbf;
    else
        u = unew;
        v = vnew;
    end

    t = kt*dt;

  if ( (mod(kt,25)==0) || (kt<=10) || (kt==numtimesteps) )
    disp([kt t]);
    sphplot = Eplot*u;
    sphplot = reshape(sphplot, size(xp));
    set(0, 'CurrentFigure', 1);
    set(Hplot, 'CData', sphplot);
    title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
    drawnow;
  end
end
t_explicit = toc;
