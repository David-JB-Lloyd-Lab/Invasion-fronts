function [F,J] = Swift_1D_per_cond(uu,nu,L,DX,ww,uT)

n = length(uu)-1;
cx= uu(n+1);
u = uu(1:n);

F = L*u + nu*u.^3 - u.^5 + cx*DX*u;
F(n+1) = ww*((DX*uT).*(u-uT));

if nargout > 1
    J = sparse(n+2);
    J(1:n,1:n) = L + spdiags(3*nu*u.^2 - 5*u.^4,0,n,n) + cx*DX;
    J(1:n,n+1) = DX*u;
    J(n+1,1:n) = (2*pi)^2*(DX*uT)'/n;
end