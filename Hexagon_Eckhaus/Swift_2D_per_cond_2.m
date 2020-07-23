function [F,J] = Swift_2D_per_cond_2(uu,p,mesh_params)

n = length(uu)-2;
cx= uu(n+1);    % phase condition parameter
cy= uu(n+2);    % phase condition parameter
u = uu(1:n);    % solution of SH equation

mu = p(1);
nu = p(2);
ky = p(3);
kx = p(4);
uT = mesh_params.uT; % template function for phase conditions

DX = kx*mesh_params.Dx;
DY = ky*mesh_params.Dy;
ww = mesh_params.ww2D;

L  = -(mesh_params.Dxx*kx^2 + mesh_params.Dyy*ky^2 + mesh_params.Ixy)^2 + mu*speye(n);

F = L*u + nu*u.^2 - u.^3 + cx*DX*u + cy*DY*u;   % SH equation
F(n+1) = ww*((DX*uT).*(u-uT));                  % phase condition
F(n+2) = ww*(((DX+DY)*uT).*(u-uT));             % phase condition

if nargout > 1 % Jacobian!
    J = sparse(n+2);
    J(1:n,1:n) = L + spdiags(2*nu*u - 3*u.^2,0,n,n) + cx*DX + cy*DY;
    J(1:n,n+1) = DX*u;
    J(n+1,1:n) = (2*pi)^2*(DX*uT)'/n;
    J(1:n,n+2) = DY*u;
    J(n+2,1:n) = (2*pi)^2*((DX+DY)*uT)'/n;
end
