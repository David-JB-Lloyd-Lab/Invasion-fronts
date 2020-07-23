function [F] = SH_2D_eig(uin,p,mesh_params)

  % Rename parameters
n = mesh_params.nx*mesh_params.ny;

mu = p(1);
nu = p(2);
ky = p(3);
kx = p(4);
  
  % Auxiliary variables
  n  = mesh_params.nx*mesh_params.ny;
  u  = uin(1:n);
  v  = uin(1+n:2*n);
  vg = 1i*uin(2*n+1:3*n);
  vgg= uin(3*n+1:4*n);
  
  cx    = uin(4*n+1);
  cy    = uin(4*n+2);
  lam  = uin(4*n+3);
  lamg = 1i*uin(4*n+4);
  lamgg= uin(4*n+5);
  
  uT = mesh_params.uT;

  DX = kx*mesh_params.Dx;
  DY = ky*mesh_params.Dy;
  ww = mesh_params.ww2D;

% set up nonlinear system to solve
  L  = -(mesh_params.Dxx*kx^2 + mesh_params.Dyy*ky^2 + mesh_params.Ixy)^2 + mu*speye(n);
  
  F1 = L*u + nu*u.^2 - u.^3 + cx*DX*u + cy*DY*u;
  
  DN = @(u) 2*nu*u - 3*u.^2;
  
  L = -(mesh_params.Dxx*kx^2 + mesh_params.Dyy*ky^2 + mesh_params.Ixy)^2 + mu*speye(n)  ...
      + spdiags(DN(u),0,n,n);
  
  F2 = L*v - lam*v;
  
  F3 = imag(L*vg - lamg*v - 4*1i*(speye(n) + kx^2*mesh_params.Dxx)*kx*mesh_params.Dx*v);

  F4 = real(L*vgg - 2*lamg*vg - lamgg*v ...
      - 8*1i*(speye(n) + kx^2*mesh_params.Dxx)*kx*mesh_params.Dx*vg ...
      + 4*v + 12*kx^2*mesh_params.Dxx*v);
  
F = [F1;F2;F3;F4];

  F(4*n+1) =  ww*((DX*uT).*(u-uT));
  F(4*n+2) =  ww*(((DX+DY)*uT).*(u-uT));
  F(4*n+3) =  norm(v,2).^2-1;
  F(4*n+4) =  ww*(v.*imag(vg));
  F(4*n+5) =  ww*(v.*vgg);
  
  

