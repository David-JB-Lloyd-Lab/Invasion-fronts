function [F,J] = Swift_2D_per_cond_nonvar(uu,nu,alpha,L,DX,DY,ww,uT)

n = length(uu)-2;
cx= uu(n+1);
cy= uu(n+2);
u = uu(1:n);

F = L*u + nu*u.^2 - u.^3 + cx*DX*u + cy*DY*u + alpha*((DX*u).^2 + (DY*u).^2);
F(n+1) = ww*((DX*uT).*(u-uT)); % phase conditions
F(n+2) = ww*(((DX+DY)*uT).*(u-uT));

if nargout > 1
    J = sparse(n+2);
    J(1:n,1:n) = L + spdiags(2*nu*u - 3*u.^2,0,n,n) + cx*DX + cy*DY ...
        +2*alpha*spdiags(DX*u,0,n,n)*DX + 2*alpha*spdiags(DY*u,0,n,n)*DY;
    J(1:n,n+1) = DX*u;
    J(n+1,1:n) = (2*pi)^2*(DX*uT)'/n;
    J(1:n,n+2) = DY*u;
    J(n+2,1:n) = (2*pi)^2*((DX+DY)*uT)'/n;
end