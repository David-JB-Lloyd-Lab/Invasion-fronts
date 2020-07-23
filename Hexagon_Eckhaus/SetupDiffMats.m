%% setup differentiation matrices
% Fourier differentiation matrix first order for x between -pi and pi
column = [0 .5*(-1).^(1:nx-1).*cot((1:nx-1)*hx/2)]';
D1x = toeplitz(column,column([1 nx:-1:2]));

% Fourier differentiation matrix for t between -pi and pi
D2x = toeplitz([-pi^2/(3*hx^2)-1/6 .5*(-1).^(2:nx)./sin(hx*(1:nx-1)/2).^2]);
D4x = D2x^2;

wx = 2*pi*ones(nx,1)/nx; % integration weights for trapzoid rule - mean

%wzx = kron(wz,wy'); wzy=wzy(:)';

Ix = speye(nx);

% Fourier differentiation matrix first order for y between -pi and pi
column = [0 .5*(-1).^(1:ny-1).*cot((1:ny-1)*hy/2)]';
D1y = toeplitz(column,column([1 ny:-1:2]));

% Fourier differentiation matrix for t between -pi and pi
D2y = toeplitz([-pi^2/(3*hy^2)-1/6 .5*(-1).^(2:ny)./sin(hy*(1:ny-1)/2).^2]);
D4y = D2y^2;

wy = 2*pi*ones(ny,1)/ny; % integration weights for trapzoid rule - mean

Iy = speye(ny); 

% diff mats

Lap = kron(D2x,Iy) + kron(Ix,D2y);
Dx  = kron(D1x,Iy); Dy = kron(Ix,D1y);
Dxx = kron(D2x,Iy); Dyy= kron(Ix,D2y);
L = -(speye(nx*ny) + Lap)^2;
Ixy = speye(nx*ny);

[ix,iy] = meshgrid(x,y);

ww2D = kron(wx,wy); ww2D = ww2D(:)';


mesh_params.nx  = nx; mesh_params.ny  = ny; 
mesh_params.Lx  = Lx; mesh_params.Ly  = Ly;

mesh_params.ix  = ix;  mesh_params.iy  = iy;
mesh_params.x   = x;   mesh_params.y   = y;

mesh_params.wx  = wx; mesh_params.wy  = wy;

mesh_params.D1y = D1y; mesh_params.D2y = D2y; mesh_params.D4y = D4y;
mesh_params.D1x = D1x; mesh_params.D2x = D2x; mesh_params.D4x = D4x;

mesh_params.Dxx = Dxx; mesh_params.Dyy = Dyy; mesh_params.Ixy = Ixy;
mesh_params.Dx  = Dx;  mesh_params.Dy  = Dy;

mesh_params.ww2D = ww2D;
