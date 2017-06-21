%if (i < 60) skip = 1; else skip = 20; end
skip = 25;

% plot every 'skip' timesteps
if (i==1) || (mod(i,skip) == 0)
  fn = 1;
  if (i==1)
    figure(fn);
  end
  set(0,'CurrentFigure', fn); clf;

  %if exist('cpx', 'var')
  %  wh = 'surf';
  if exist('z', 'var')
    wh = '3d';
  elseif exist('yy', 'var')
    wh = '2d';
  else
    wh = '1d';
  end

  switch wh
    case '1d'
      plot(x, u, 'kx-', 'linewidth', 3);
      hold on
      plot(x, v, 'r-', 'linewidth', 3);
      xlabel('x'); ylabel('u');
      legend('u', 'v')
    case '2d'
      up = reshape(u, N, N);
      pcolor(xx,yy,up);
      axis equal; axis tight
      shading flat
      %colorbar
      xlabel('x'); ylabel('y');
    case '3d'
      up = reshape(u, N, N, N);
      isosurface(x,y,z,up,0.45);
      axis equal; %axis tight
      axis([x(1) x(end) y(1) y(end) z(1) z(end)]);
    case 'surf'
      if i==1
        [th, phi] = meshgrid(0:pi/64:2*pi, -pi/2:pi/64:pi/2);
        [xp,yp,zp] = sph2cart(th, phi, Radius);
      end
      %up = ba_interp3(x,y,z, reshape(u,N,N,N), xp,yp,zp, 'cubic');
      % ba_interp3 is faster, but below works too:
      up = interp3(x,y,z, reshape(u,N,N,N), xp,yp,zp);
      surf(xp, yp, zp, up)
      %colorbar
      axis equal
      shading flat
      camlight left
  end
  t = i*dt;
  title(['u, t = ', num2str(t)]);
  drawnow()

end
