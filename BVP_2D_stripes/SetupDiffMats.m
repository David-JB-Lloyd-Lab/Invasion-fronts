%% setup differentiation matrices

% Differentiation matrix in z: finite differences
ez = ones(nz,1); 
Dz = spdiags([ez -8*ez 0*ez 8*ez -ez],-2:2, nz, nz);
Dz(1,:)   = 0; Dz(2,2)   = 1;
Dz(nz,:) = 0; Dz(nz-1,nz-1) = -1;
Dz = Dz/(12*hz);

D2z = spdiags([-ez 16*ez -30.*ez 16*ez -ez], -2:2, nz, nz);
D2z(1,2)= 32; D2z(1,3)= -2;
D2z(2,1)= 16; D2z(2,2)= -31; 
D2z(nz,:) = D2z(1,nz:-1:1);
D2z(nz-1,:) = D2z(2,nz:-1:1);
D2z = D2z/(12*hz^2);

D4z = spdiags([-ez 12*ez -39*ez 56*ez -39*ez +12*ez -ez],-3:3, nz, nz);
D4z(1,2) = -78; D4z(1,3) = 24; D4z(1,4) = -2;
D4z(2,2) =  68; D4z(2,3) = -40;
D4z(3,2) = -40;
D4z(nz,:) = D4z(1,nz:-1:1);
D4z(nz-1,:) = D4z(2,nz:-1:1);
D4z(nz-2,:) = D4z(3,nz:-1:1);
D4z = D4z/(6*hz^4);

Iz = speye(nz);

% wz = hx*[1, 2*ones(1,nz-2)+2*mod([1:nz-2],2),1]/3;   % Simpson weights for intergration int = w*u un-scaled  must *hx
wz = ones(nz,1)/nz;

% Fourier differentiation matrix first order for t between -pi and pi
column = [0 .5*(-1).^(1:nt-1).*cot((1:nt-1)*ht/2)]';
D1t = toeplitz(column,column([1 nt:-1:2]));

% Fourier differentiation matrix for t between -pi and pi
D2t = toeplitz([-pi^2/(3*ht^2)-1/6 .5*(-1).^(2:nt)./sin(ht*(1:nt-1)/2).^2]);
D4t = D2t^2;

wt = 2*pi*ones(nt,1)/nt; % integration weights for trapzoid rule - mean

wzt = kron(wz,wt'); wzt=wzt(:)';

It = speye(nt); 

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

wzy = kron(wz,wy'); wzy=wzy(:)';

Iy = speye(ny); 

% diff mats
DZ = sparse(kron(It, kron(Dz,  Iy)));
DZZ= sparse(kron(It, kron(D2z, Iy)));
D4Z= sparse(kron(It, kron(D4z, Iy)));
DT = sparse(kron(Iy, kron(Iz, D1t)));
DTT= sparse(kron(Iy, kron(Iz, D2t)));
DY = sparse(kron(D1y,kron(Iz, It )));
DYY= sparse(kron(D2y,kron(Iz, It )));
D4Y= sparse(kron(D4y,kron(Iz, It )));
D2Z2Y= DZZ*DYY;

Lap = kron(D2x,Iy) + kron(Ix,D2y);
Dx  = kron(D1x,Iy); Dy = kron(Ix,D1y);
Dxx = kron(D2x,Iy); Dyy= kron(Ix,D2y);
L = -(speye(nx*ny) + Lap)^2;
Ixy = speye(nx*ny);

[ix,iy] = meshgrid(x,y);

ww   = kron(kron(wz,wt),wy); ww = ww(:)';
ww2D = kron(wx,wy); ww2D = ww2D(:)';

[zz,tt,yy] = meshgrid(z,t,y);
[iz,it] = meshgrid(z,t);

mesh_params.nz  = nz;  mesh_params.nt  = nt; mesh_params.nx  = nx; mesh_params.ny  = ny; 
mesh_params.Lz  = Lz;  mesh_params.Lt  = Lt; mesh_params.Lx  = Lx; mesh_params.Ly  = Ly;

mesh_params.zz  = zz;  mesh_params.tt  = tt; mesh_params.yy  = yy;
mesh_params.ix  = ix;  mesh_params.iy  = iy;
mesh_params.iz  = iz;  mesh_params.it  = it;
mesh_params.x   = x;   mesh_params.y   = y;

mesh_params.Iz  = Iz;  mesh_params.It  = It; mesh_params.Iy  = Iy; mesh_params.Ix = Ix;
mesh_params.wz  = wz;  mesh_params.wt  = wt; mesh_params.wy  = wy;

mesh_params.Dz  = Dz; mesh_params.D2z = D2z; mesh_params.D4z = D4z;
mesh_params.D1t = D1t; mesh_params.D2t = D2t; mesh_params.D4t = D4t;
mesh_params.D1y = D1y; mesh_params.D2y = D2y; mesh_params.D4y = D4y;
mesh_params.D1x = D1x; mesh_params.D2x = D2x; mesh_params.D4x = D4x;

mesh_params.Dxx = Dxx; mesh_params.Dyy = Dyy; mesh_params.Ixy = Ixy;
mesh_params.Dx  = Dx;  mesh_params.Dy  = Dy;

mesh_params.DZZ = DZZ; mesh_params.DTT = DTT; mesh_params.DYY = DYY; 
mesh_params.DZ = DZ;   mesh_params.DT = DT; mesh_params.DY = DY;
mesh_params.D2Z2T = D2Z2Y;
mesh_params.D4Z = D4Z; mesh_params.D4Y = D4Y;
mesh_params.ww = ww; mesh_params.ww2D = ww2D;mesh_params.wx = wx';
mesh_params.iend=find(abs(zz(:))>Lz/2-4*hx);

