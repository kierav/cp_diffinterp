clear all
%% Heat equation on a circle
% This example solves the heat equation on a 2D circle, with initial
% conditions u = cos(theta), and exact solution u(t) =
% exp(-t)*cos(theta)


% adjust as appropriate
%addpath('../cp_matrices');
%addpath('../surfaces');

useLocal = 1; % 1 - use 4x4 points, 0 - use all points

%% Construct a grid in the embedding spaces
ep=1;
dx = 0.1;                   % grid size

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;

nx = length(x1d);
ny = length(y1d);

%% Find closest points on the surface
% For each point (x,y), we store the closest point on the circle
% (cpx,cpy)

% meshgrid is only needed for finding the closest points, not afterwards
[xx yy] = meshgrid(x1d, y1d);
% function cpCircle for finding the closest points on a circle
[cpx, cpy, dist] = cpCircle(xx,yy);
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


%% Function u in the embedding space
% u is a function defined on the grid (eg heat if solving the heat
% equation)
[th, r] = cart2pol(xg,yg);
u = cos(th);
initialu = u;       % store initial value


%% Construct an interpolation matrix for closest point
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');

E = interp2_matrix(x1d, y1d, cpxg, cpyg, p, band);
[Ei,Ej,Es] = interp2_matrix(x1d,y1d,cpxg,cpyg,p,band);
Ej = reshape(Ej,length(cpxg),(p+1)^2);
% e.g., closest point extension:
%u = E*u;

%% Create Laplacian matrix for heat equation

L = laplacian_2d_matrix(x1d,y1d, order, band);


%% Construct an interpolation matrix for plotting on circle

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,100)';
r = ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix(x1d, y1d, xp, yp, p, band);


figure(1); clf;
figure(2); clf;
figure(3); clf;


%% Time-stepping for the heat equation
Tf = 2;
dt = 0.2*dx^2;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;

% construct rbf matrix
D = zeros(length(cpxg));
if useLocal == 1
    for j = 1:length(cpxg)
        x = xg(Ej(j,:));
        y = yg(Ej(j,:));
        [A,B] = rbf(ep,cpxg(j),cpyg(j),x,y);
    %     D(j,Ej(j,:)) = B*pinv(A);
        D(j,Ej(j,:)) = B/A;
    %     [A,B] = rbf(1,cpxg(j),cpyg(j),xg,yg);
    %     D(j,:) = B*pinv(A);
    end
else

    [A,B] = rbf(ep,cpxg,cpyg,xg,yg);
    D = B*pinv(A);   
end

% D(1,:) = 0;
% D(end,:) = 0;
% D(1,end) = 1;
% D(end,1) = 1;


%%plots the stability region of FE
%figure
%plot(real(eig(D)),imag(eig(D)),'*')
%pause 
%hold on
%z = exp(1i*pi*(0:200)/100); r = z-1;
%plot(r/dt, 'r')
%pause
   
% Trying to check our D matrix
uexactdiff = @(theta) -1*cos(theta);
[th1, r1] = cart2pol(cpxg,cpyg);
ucheck = uexactdiff(th1);
plot(th1,D*initialu,'x');
%plot3(cpxg,cpyg, D*initialu, 'x');
hold on;
plot(th1,ucheck,'rx');
%plot3(cpxg,cpyg, ucheck, 'rx');
pause

for kt = 1:numtimesteps
  % explicit Euler timestepping
  unew = u + dt*D*u;
%   u = unew;
  % closest point extension
  if ( mod(kt,10) == 0)
     u = E*unew;
  else
     u = unew;
  end

  t = kt*dt;

  if ( (kt < 10) || (mod(kt,10) == 0) || (kt == numtimesteps) )
    %% plot over computation band
    %plot2d_compdomain(u, xg, yg, dx, dx, 1)
    %title( ['embedded domain: soln at time ' num2str(t) ...
    %        ', timestep #' num2str(kt)] );
    %xlabel('x'); ylabel('y');
    %%hold on
    %plot(xp,yp,'k-', 'linewidth', 2);
    %%axis equal;  axis tight;
    %circplot = Eplot*u;
    set(0, 'CurrentFigure', 1);
    clf;
    exactplot = exp(-t)*cos(th1);
    plot3(cpxg,cpyg, u, 'x');
    title( ['soln at time ' num2str(t) ', on circle'] );
    %xlabel('theta'); ylabel('u');
    hold on;
    % plot analytic result
    plot3(cpxg,cpyg, exactplot, 'ro');
    
    % plot value on circle
    set(0, 'CurrentFigure', 2);
    clf;
    circplot = Eplot*u;
    exactplot = exp(-t)*cos(thetas);
    plot(thetas, circplot);
    title( ['soln at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('u');
    hold on;
    % plot analytic result
    plot(thetas, exactplot, 'r--');
    legend('explicit Euler', 'exact answer', 'Location', 'SouthEast');
    error_circ_inf = max(abs( exactplot - circplot ));

    set(0, 'CurrentFigure', 3);
    clf;
    plot(thetas, circplot - exactplot);
    title( ['error at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('error');

    %pause
    drawnow();
  end
end
