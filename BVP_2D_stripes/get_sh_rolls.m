function [uu_r,Duu_r,ur,ky_new] = get_sh_rolls(p,kz,ky,mesh_params,ur_0)

  mu = p(1); 
  nu = p(2); 
  % construct Swift-Hohenberg operators
  Dx = mesh_params.D1x*kz;
  L  = -(mesh_params.D2x*kz^2 + mesh_params.Ix)^2 + mu*mesh_params.Ix;
  ww = mesh_params.wx;
  
  % Solve SH on a periodic domain find rolls
   sh_rhs_1D = @(u) Swift_1D_per_cond(u,nu,L,Dx,ww,mesh_params.uT);
   options = optimset('Jacobian','on','Display','off');
   [ur] = fsolve(sh_rhs_1D,[ur_0;0],options);
   if (abs(ur(end))>1e-3)
       disp(['Warning phase condition for rolls is non-zero =',num2str(ur(end-1))]);
   end

  ur = ur(1:end-1); % remove the wave speed from ur
  Dur = Dx*ur;
  
  uur = repmat(ur,1,mesh_params.ny);
  Duur= repmat(Dur,1,mesh_params.ny);
  
  % Interpolating on the eta mesh
  eta = kz*mesh_params.iz(:) + mesh_params.it(:);  % find skewed rolls and interpolate
  xi  = mesh_params.y(:);
  
  uu_r = zeros(size(eta));
  Duu_r= zeros(size(eta));
  
  uu_r = periodic_interp_2D(uur,mesh_params.x,mesh_params.y,eta,xi);
  Duu_r= periodic_interp_2D(Duur,mesh_params.x,mesh_params.y,eta,xi);

  uu_r = uu_r(:);
  Duu_r= Duu_r(:);

end

function p = periodic_interp_2D(f,xn,yn,xf,yf)

[M,N]=size(f);

xn=xn(:); xf=xf(:); yn=yn(:); yf=yf(:);
Mf=length(xf);      Nf=length(yf);

if or(length(xn)~=M,length(yn)~=N)
    error('grid and data size do not match');
end

% Get distances between nodes and interpolation points
xdist=repmat(xf,1,M)-repmat(xn.',Mf,1);
ydist=repmat(yf,1,N)-repmat(yn.',Nf,1);

Hx = diric(xdist,M).*cos(xdist/2);
Hy = diric(ydist,N).*cos(ydist/2);
p  = Hx*f*Hy.';

end
