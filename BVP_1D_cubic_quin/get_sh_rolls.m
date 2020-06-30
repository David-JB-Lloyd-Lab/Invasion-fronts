function [uu_r,Duu_r,ur] = get_sh_rolls(p,kz,mesh_params,ur_0)

  mu = p(1); 
 k = kz;
  
  % construct Swift-Hohenberg operators
  Dt_scale = mesh_params.Dt*k;
  L_1D     = -(mesh_params.D4t*k^4 + 2*mesh_params.D2t*k^2 + mesh_params.It);
  
  % Solve SH on a periodic domain 
  nt = mesh_params.nt; ht = 2*pi/nt;  t = ht*(1:nt);
  
  % find 1D rolls
  sh_rhs_1D = @(u) rhs_sh_phaseCond(u,p,L_1D,Dt_scale,mesh_params.tt(:,1),ur_0);
  options = optimset('Jacobian','on','Algorithm','levenberg-marquardt','Display','off');
  [ur] = fsolve(sh_rhs_1D,[ur_0;0],options);
  
  if abs(ur(end))>1e-3
      disp(['Warning phase condition for rolls is non-zero =',num2str(ur(end))]);
  end
  
  ur = ur(1:end-1); % remove the translation speed from ur (this should be zero!)
  
  Dur = Dt_scale*ur;
  
  % Interpolating on the eta mesh
  eta = kz*mesh_params.zz(:) + mesh_params.tt(:);  % find skewed rolls and interpolate

  uu_r = zeros(size(eta));
  Duu_r= zeros(size(eta));
  periodic_sinc = ones(size(eta)); 

for i = 1:nt
    zz = eta-t(i);
    periodic_sinc = diric(zz,nt).*cos(zz/2);
    uu_r = uu_r + ur(i)*periodic_sinc;
    Duu_r= Duu_r+Dur(i)*periodic_sinc;
end

