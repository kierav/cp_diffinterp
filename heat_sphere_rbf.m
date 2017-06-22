%% Heat equation on a sphere
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.
%
% This example solves the heat equation on the surface of a sphere,
% with initial conditions u = cos(4*theta)

%%
useLocal = 1; % 1 - use 4x4x4 points, 0 - use all points

% 3D example on a sphere
% Construct a grid in the embedding space
Er_final=[];
for ep =.66:0.01:0.69%.65:0.05:1%0.05:0.05:1;
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
[cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);
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


%% Function u in the embedding space
% u is a function defined on the grid (eg heat if solving the heat
% equation)

% assign some initial value (using initial value of cos (8*theta))
[th, phi, r] = cart2sph(xx,yy,zz);
u = cos(phi + pi/2);

% this makes u into a vector, containing only points in the band
u = u(band);

initialu = u;       % store initial value


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

% D(1,:) = 0;
% D(end,:) = 0;
% D(1,end) = 1;
% D(end,1) = 1;



   
% % Trying to check our D matrix
% uexactdiff = @(theta) -1*cos(theta);
% [th1, r1] = cart2pol(cpxg,cpyg);
% ucheck = uexactdiff(th1);
% plot(th1,D*initialu,'x');
% %plot3(cpxg,cpyg, D*initialu, 'x');
% hold on;
% plot(th1,ucheck,'rx');
% %plot3(cpxg,cpyg, ucheck, 'rx');
% pause

%% Construct an interpolation matrix for plotting on sphere

% plotting grid on sphere, based on parameterization
[xp,yp,zp] = paramSphere(64);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, r] = cart2sph(xp1,yp1,zp1);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);



%% Time-stepping for the heat equation

Tf = 1;
dt = 1/6*dx^2;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;
%plots the stability region of FE
figure(2);clf;
plot(real(eigs(D)),imag(eigs(D)),'*')
%pause 
hold on
z = exp(1i*pi*(0:200)/100); r = z-1;
plot(r/dt, 'r')
%pause
tic
for kt = 1:numtimesteps
    % explicit Euler timestepping
    unew = u + dt*(D*u);

    % closest point extension
    if ( mod(kt,5) == 0)
        u = E*unew;
    else
        u = unew;
    end

    t = kt*dt;

    % plot value on sphere
    if (mod(kt,100) == 0) || (kt < 10) || (kt == numtimesteps)
      %figure(1);
      sphplot = Eplot*u;

	  err = norm(exp(-2*t)*cos(phi_plot + pi/2)-sphplot,inf) / norm(exp(-2*t)*cos(phi_plot + pi/2),inf);
      [t dt dx err]

      sphplot = reshape(sphplot, size(xp));
      figure(1); clf;
      surf(xp, yp, zp, sphplot);
      title( ['soln at time ' num2str(t) ', kt= ' num2str(kt)] );
      xlabel('x'); ylabel('y'); zlabel('z');
      %caxis([-1.05 1.05]);   % lock color axis
      axis equal; shading interp;
%      if ~exist OCTAVE_VERSION camlight left;
      colormap(jet);
      colorbar;
      drawnow();% pause(0);
%      end
    end
end
t_explicit = toc;
Er_final=[Er_final,err];
end
figure
semilogy(ep,Er_final,'.')
title('Using IMQ for heat equ on Sphere')
xlabel('\epsilon')
ylabel('error')