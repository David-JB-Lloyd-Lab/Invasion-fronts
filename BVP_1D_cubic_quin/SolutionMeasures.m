function F = SolutionMeasures(step,u,p,mesh_params)
  
global ur0_left;

  % Solution measures (they are displayed on screen)
  % by setting stepperPars.PlotBranchVariableId = k, 
  % the kth solution measure is plotted in the 
  % bifurcation diagram on the fly
  
       % Rename parameters
  mu    = p(1);
  nu    = p(2);
  m     = p(5);
  d     = p(6);
  
  % Auxiliary variables
  n = mesh_params.nz*mesh_params.nt;
  w     = u(1:n);
  kzm   = u(n+1); 
  kt    = u(n+2); 
 
  chi_p = 1/2 + 1/2*tanh(m*(mesh_params.zz(:)-d));
  chi_m = chi_p(end:-1:1);
  Dchi_p = mesh_params.DZ * chi_p;
  Dchi_m = mesh_params.DZ * chi_m;

  % find skewed rolls
  [uu_rm,Duu_rm,~] = get_sh_rolls(p,kzm,mesh_params,ur0_left);
   
   nz = mesh_params.nz;
   nt = mesh_params.nt;
   w_out = u(1:nz*nt) + chi_m.*uu_rm;
   
   u2 = reshape(u(1:nz*nt),nt,nz);
 
      % detection of pitchfork bifurcation -> compute det(J)
    [F,J]= Front_invasion_SH_phase_condition(u,p,mesh_params,0);
    B = sign(det(J(1:end-2,1:end-2)));

    D1 = mesh_params.wzt*(mesh_params.DZ*w_out).^2 - mesh_params.wzt*(mesh_params.D2Z*w_out).^2;
   
    [V,D] = eigs(J(1:end-2,1:end-2)',2,0.1);
    [~,ii] = sort(real(diag(D)));
    phi = V(:,ii(1));
    uz  = mesh_params.DZ*w_out;
    uzzz= mesh_params.DZ*mesh_params.D2Z*w_out;
    trans_lam1 = 2*mesh_params.wzt*((uz + uzzz).*phi)./(mesh_params.wzt*(kt*uz.*phi));
    
    phi = V(:,ii(2));
    uz  = mesh_params.DZ*w_out;
    uzzz= mesh_params.DZ*mesh_params.D2Z*w_out;
    trans_lam2 = 2*mesh_params.wzt*((uz + uzzz).*phi)./(mesh_params.wzt*(kt*uz.*phi));
    
    dd = eigs(J(1:end-2,1:end-2)',5,0.1);
    
   F = [kzm kt B D1 trans_lam1 trans_lam2];
  
  

end
