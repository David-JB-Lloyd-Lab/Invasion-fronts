function [F,J] = rhs_sh_phaseCond(uu,p,L,Dt,y,uT)

  % Rename parameters
  mu = p(1); 
  nu = p(2);
  n = size(uu,1)-1;
  u = uu(1:n);
  c = uu(n+1);
  

  % Right-hand side
%   uT = temp;
  F = zeros(n+1,1);
  F(1:n) =  L*u + mu*u + nu*u.^3 - u.^5 + c*Dt*u;
  F(n+1) = mean( (Dt*uT) .* ( u - uT)) ;

  % Jacobian
  if nargout > 1
     J = sparse(n+1);
     J(1:n,1:n) = L + spdiags( +mu + 3*nu*u.^2 - 5*u.^4, 0, n, n) + c*Dt;
     J(:,n+1)   = Dt*u;
     J(n+1,1:n) = (Dt*uT)'/n;
  end

end
