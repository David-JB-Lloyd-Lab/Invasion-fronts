function [F,J] = Front_invasion_SH_phase_condition_nonvar(u,p,mesh_params,jac)

  %% Setup
  global ur0_left;

  % Rename parameters
  mu    = p(1);
  nu    = p(2);
  m     = p(5);
  d     = p(6);
  alpha = p(8);
  
  % Auxiliary variables
  n = mesh_params.nz*mesh_params.nt*mesh_params.ny;
  w     = u(1:n); % core function
  kzm   = u(n+1); % far-field selected wavenumber kzm
  kt    = u(n+2); % selected invasion wavenumber
  ky = p(7);
  c     = kt/kzm; % front wave speed
  
  chi_m = 1/2 + 1/2*tanh(m*(-mesh_params.zz(:)-d));

  % find skewed rolls
  if jac ~= 1
    [uu_rm,Duu_rm,ur0_left ] = get_sh_hexs_nonvar(p,kzm,ky,mesh_params,ur0_left);
  else
    [uu_rm,Duu_rm,~]         = get_sh_hexs_nonvar(p,kzm,ky,mesh_params,ur0_left);
  end
  
  %% Right-hand side
  % Linear operators
  L = -mesh_params.D4Z - ky^4*mesh_params.D4Y ...
      - 2*ky^2*mesh_params.DZZ*mesh_params.DYY ...
      - 2*mesh_params.DZZ - 2*ky^2*mesh_params.DYY ... 
      - spdiags((1-mu)*ones(n,1),0,n,n) ...
      +c*mesh_params.DZ - kt*mesh_params.DT; 
  
  L2 = -mesh_params.D4Z - ky^4*mesh_params.D4Y ...
      - 2*ky^2*mesh_params.DZZ*mesh_params.DYY ...
      - 2*mesh_params.DZZ - 2*ky^2*mesh_params.DYY ... 
      - spdiags((1-mu)*ones(n,1),0,n,n);
  
  CHI_m = spdiags(chi_m,0,n,n);
  
  % nonlinearity
  N  = @(u) nu*u.^2 - u.^3 + alpha*((mesh_params.DZ*u).^2 + (ky*mesh_params.DY*u).^2);
  
  F = zeros(length(u),1);
  
  % main core PDE problem
  F(1:n) = L*(uu_rm.*chi_m + w ) + N(uu_rm.*chi_m + w ) ...
      - chi_m.*(L2*uu_rm + N(uu_rm)); 
  
  F(mesh_params.iend) = w(mesh_params.iend); % dirchlet bcs at the ends of the domain
  
  % phase condition for phase field rolls
  loc_wold = mesh_params.loc;
  F(n+1) = mesh_params.ww*(loc_wold.*(w  - mesh_params.w0));

  % Phase condition for front speed int w^{old}_z(w-w^{old})
  iBm = find( mesh_params.zz(:) <= -mesh_params.Lz/2 + 2*pi/abs(kzm) ); % find indices close to x = Lx - 2*pi/kx : Lx and y = 0..Ly
  z = mesh_params.z; hz= z(2) - z(1);
  iim= find(z <= -mesh_params.Lz/2 + 2*pi/abs(kzm)); nnzm = length(iim);
  
  % phase condition 2: find du_r
  uuu_rm_prime = Duu_rm;
  uuu_rm_prime = uuu_rm_prime(iBm); % define uuu_r_prime for integral
  
  wz = [1, 2*ones(1,nnzm-2)+2*mod([1:nnzm-2],2),1]*hz/3;   % Simpson weights for intergration int = w*u
  wt = 2*mesh_params.Lt*ones(mesh_params.nt,1)/mesh_params.nt;
  wy = 2*mesh_params.Ly*ones(mesh_params.ny,1)/mesh_params.ny;
  wwm= kron(kron(wz,wt),wy);
  wwm= wwm(:)';
  
  wBm = w(iBm);                       % find w on the domain x = Lx - 2*pi/kx : Lx and y = 0..Ly
  F(n+2) = wwm*(uuu_rm_prime .* wBm);  % phase conditon \int w* du_r/dkx dx

  %% Jacobian
  if nargout > 1
     J = sparse(n,n);
     
     DN = @(u) spdiags(2*nu*u - 3*u.^2,0,n,n) ...
         + 2*alpha*(spdiags(mesh_params.DZ*u,0,n,n)*mesh_params.DZ + ky*spdiags(ky*mesh_params.DY*u,0,n,n)*mesh_params.DY);
     
     J(1:n,1:n) = L + DN(uu_rm.*chi_m + w);
      
     J(mesh_params.iend,:) = 0;
     J(mesh_params.iend,mesh_params.iend)=speye(length(mesh_params.iend));

     J(n+1,1:n)  =mesh_params.ww*spdiags(loc_wold,0,n,n);
     J(n+2,iBm)  = wwm*spdiags(uuu_rm_prime,0,nnzm*mesh_params.nt*mesh_params.ny,nnzm*mesh_params.nt*mesh_params.ny);

     epsiF = 1e-8;
     dF = Front_invasion_SH_phase_condition_nonvar([u(1:n); u(n+1) + epsiF; u(n+2)],p,mesh_params,1);
     J(:,n+1) = (dF - F)/epsiF;

     epsiF = 1e-8;
     dF = Front_invasion_SH_phase_condition_nonvar([u(1:n); u(n+1); u(n+2) + epsiF],p,mesh_params,1);
     J(:,n+2) = (dF - F)/epsiF;
     
  end
      
end
