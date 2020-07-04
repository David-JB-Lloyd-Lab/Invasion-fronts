function [F,J] = Front_invasion_SH_phase_condition(u,p,mesh_params,jac)

global ur0_left;

  % Rename parameters
  mu    = p(1);
  nu    = p(2);
  m     = p(5);
  d     = p(6);
  
  % Auxiliary variables
  n = mesh_params.nz*mesh_params.nt*mesh_params.ny;
  w     = u(1:n);
  kzm   = u(n+1); 
  kt    = u(n+2); 
  cy    = u(n+3);
  ky = p(7);
  c     = kt/kzm; % front wave speed
  %
  chi_p = 1/2 + 1/2*tanh(m*(mesh_params.zz(:)-d));
  chi_m = chi_p(end:-1:1);
  Dchi_m = mesh_params.DZ * chi_m;

  % find skewed rolls
  if jac ~= 1
    [uu_rm,Duu_rm,ur0_left ] = get_sh_rolls(p,kzm,ky,mesh_params,ur0_left);
  else
    [uu_rm,Duu_rm,~]         = get_sh_rolls(p,kzm,ky,mesh_params,ur0_left);
  end
  % Right-hand side
  L = -mesh_params.D4Z - ky^4*mesh_params.D4Y ...
      - 2*ky^2*mesh_params.DZZ*mesh_params.DYY ...
      - 2*mesh_params.DZZ - 2*ky^2*mesh_params.DYY ... 
      - spdiags((1-mu)*ones(n,1),0,n,n) ...
      +c*mesh_params.DZ - kt*mesh_params.DT ...
      +cy*ky*mesh_params.DY; 
  
  L2 = -mesh_params.D4Z - ky^4*mesh_params.D4Y ...
      - 2*ky^2*mesh_params.DZZ*mesh_params.DYY ...
      - 2*mesh_params.DZZ - 2*ky^2*mesh_params.DYY ... 
      - spdiags((1-mu)*ones(n,1),0,n,n);
  
  CHI_m = spdiags(chi_m,0,n,n);
  
  N  = @(u) nu*u.^3 - u.^5 ;
  
  F = zeros(length(u),1);
  
  F(1:n) = L*(uu_rm.*chi_m + w ) + N(uu_rm.*chi_m + w ) ...
      - chi_m.*(L2*uu_rm + N(uu_rm));
  
  F(mesh_params.iend) = w(mesh_params.iend); % dirchlet bcs at the ends of the domain
  
  % phase condition for phase field rolls
  loc_wold = mesh_params.loc_wold;
  F(n+1) = mesh_params.ww*(loc_wold.*(w  - mesh_params.w0));

  % Phase condition for front speed int w^{old}_z(w-w^{old})
  iBm = find( mesh_params.zz(:) <= -mesh_params.Lz/2 + 2*pi/abs(kzm) ); % find indices close to x = Lx - 2*pi/kx : Lx and y = 0..Ly

  z = linspace(-mesh_params.Lz/2,mesh_params.Lz/2,mesh_params.nz);  hz= z(2) - z(1);
  iim= find(z <= -mesh_params.Lz/2 + 2*pi/abs(kzm)); nnzm = length(iim);
  
  % phase condition 2: find du_r
  uuu_rm_prime = Duu_rm;

  % define uuu_r_prime for integral
  uuu_rm_prime = uuu_rm_prime(iBm);
  
  wz = [1, 2*ones(1,nnzm-2)+2*mod([1:nnzm-2],2),1]*hz/3;   % Simpson weights for intergration int = w*u
  wt = 2*mesh_params.Lt*ones(mesh_params.nt,1)/mesh_params.nt;
  wy = 2*mesh_params.Ly*ones(mesh_params.ny,1)/mesh_params.ny;
  wwm= kron(kron(wz,wt),wy);
  wwm= wwm(:)';
  
  wBm = w(iBm);                       % find w on the domain x = Lx - 2*pi/kx : Lx and y = 0..Ly
  F(n+2) = wwm*(uuu_rm_prime .* wBm);  % phase conditon \int w* du_r/dkx dx

    % phase condition for phase field rolls
 loc_wold_2 = (mesh_params.DY)*mesh_params.w0;
  F(n+3) = mesh_params.ww*(loc_wold_2.*(w  - mesh_params.w0));
  
  % Jacobian
  if nargout > 1
     J = sparse(n,n);
     
     DN = @(u) 3*nu*u.^2 - 5*u.^4 ;
     
     J(1:n,1:n) = L + spdiags(DN(uu_rm.*chi_m + w), 0, n, n);
      
     J(mesh_params.iend,:) = 0;
     J(mesh_params.iend,mesh_params.iend)=speye(length(mesh_params.iend));

     J(n+1,1:n)  =mesh_params.ww*spdiags(loc_wold,0,n,n);
     J(n+2,iBm)  = wwm*spdiags(uuu_rm_prime,0,nnzm*mesh_params.nt*mesh_params.ny,nnzm*mesh_params.nt*mesh_params.ny);
     J(n+3,1:n)  =mesh_params.ww*spdiags(loc_wold_2,0,n,n);

     epsiF = 1e-8;
     dF = Front_invasion_SH_phase_condition([u(1:n); u(n+1) + epsiF; u(n+2); u(n+3)],p,mesh_params,1);
     J(:,n+1) = (dF - F)/epsiF;

     epsiF = 1e-8;
     dF = Front_invasion_SH_phase_condition([u(1:n); u(n+1); u(n+2) + epsiF; u(n+3)],p,mesh_params,1);
     J(:,n+2) = (dF - F)/epsiF;
     
     epsiF = 1e-8;
     dF = Front_invasion_SH_phase_condition([u(1:n); u(n+1); u(n+2); u(n+3) + epsiF],p,mesh_params,1);
     J(:,n+3) = (dF - F)/epsiF;
  end
      
end
