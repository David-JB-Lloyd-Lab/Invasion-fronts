function F = SolutionMeasures(step,u,p,mesh_params)
  
global ur0_left;
  % Rename variables
%   n = size(u,1);
%   kxm = u(n-1);
%   kxp = u(n);
%  

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
  n = mesh_params.nz*mesh_params.nt*mesh_params.ny;
  w     = u(1:n);
  kzm   = u(n+1); 
  kt    = u(n+2);
  ky    = u(n+3);
 
%   chi_p = 1/2 + 1/2*tanh(m*(mesh_params.zz(:)-d));
%   chi_m = chi_p(end:-1:1);
%   Dchi_p = mesh_params.DZ * chi_p;
%   Dchi_m = mesh_params.DZ * chi_m;
% 
%   % find skewed rolls
%   [uu_rm,Duu_rm,~] = get_sh_rolls_2(p,kzm,ky,mesh_params,ur0_left);
%    
% %    nz = mesh_params.nz;
% %    nt = mesh_params.nt;
%    w_out = u(1:n) + chi_m.*uu_rm;
%    
%    u2 = reshape(u(1:n),mesh_params.nt,mesh_params.nz,mesh_params.ny);
 
   
% Energy 
%    E_zigzag = -1.6724205271E-01; % computed from AUTO   
%    Lap = mesh_params.D2X + mesh_params.D2Y*ky^2;
% 
%    E = mesh_params.wxy*(0.5*(w_out + Lap*w_out).^2 - 0.5*mu*w_out.^2 + w_out.^4/4 - E_zigzag);
%    
% pitchfork measure
%    u2_s = circshift(u2,ny/2,1); % compute w(x,y+pi)
%    u2_s = u2_s(:);
%    P = 0.5*mesh_params.wxy*(u(1:nx*ny) + u2_s);

   % detection of pitchfork bifurcation -> compute det(J)
%     [F,J]= Front_invasion_SH_phase_condition_2(u,p,mesh_params,0);
%     B = sign(det(J(1:end-3,1:end-3)));
%     B =  sign(det(J));
%     B = sign(det(J(1:end-2,1:end-2))-speye(nz*nt));
% D1 = mesh_params.wzt*(mesh_params.DZ*w_out).^2 - mesh_params.wzt*(mesh_params.D2Z*w_out).^2;
   
%     [V,D] = eigs(J(1:end-2,1:end-2)',2,0.1);
%     [~,ii] = sort(real(diag(D)));
%     phi = V(:,ii(1));
%     uz  = mesh_params.DZ*w_out;
%     uzzz= mesh_params.DZ*mesh_params.D2Z*w_out;
%     trans_lam1 = 2*mesh_params.wzt*((uz + uzzz).*phi)./(mesh_params.wzt*(kt*uz.*phi));
%     
%     phi = V(:,ii(2));
%     uz  = mesh_params.DZ*w_out;
%     uzzz= mesh_params.DZ*mesh_params.D2Z*w_out;
%     trans_lam2 = 2*mesh_params.wzt*((uz + uzzz).*phi)./(mesh_params.wzt*(kt*uz.*phi));
    
%     dd = eigs(J(1:end-3,1:end-3),5,0.1);
    

%   F = [mesh_params.wxy*(u(1:n-2).^2)/mesh_params.Lx abs(kxp)-abs(kxm) kxm kxp];
%   F = [mesh_params.wxy*(w_out.^2)/mesh_params.Lx abs(kxp)-abs(kxm) kxm kxp];
   F = [kzm kt ky];
  
  

end
