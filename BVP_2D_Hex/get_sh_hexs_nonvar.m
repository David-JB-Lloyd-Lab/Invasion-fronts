function [uu_r,Duu_r,ur,Dyuu_r] = get_sh_hexs_nonvar(p,kz,ky,mesh_params,ur_0)

  mu = p(1); 
  nu = p(2); 
  alpha = p(8);
  % construct Swift-Hohenberg operators
  Dx = mesh_params.Dx*kz;
  Dy = mesh_params.Dy*ky;
  L  = -(mesh_params.Dxx*kz^2 + mesh_params.Dyy*ky^2 + mesh_params.Ixy)^2 + mu;
  ww = mesh_params.ww2D;
  
  % Solve SH on a periodic domain find hexagons with ky fixed
   sh_rhs_2D = @(u) Swift_2D_per_cond_nonvar(u,nu,alpha,L,Dx,Dy,ww,mesh_params.uT);
   options = optimset('Jacobian','on','Display','off');
   [ur] = fsolve(sh_rhs_2D,[ur_0;0;0],options);
   if or((abs(ur(end))>1e-3),(abs(ur(end-1))>1e-3))
       disp(['Warning phase condition for hexagons is non-zero =',num2str(ur(end-1))]);
   end

  ur = ur(1:end-2); % remove the wave speed from ur
  Dur = Dx*ur;      % calculate derivative of solution in x
  Dyur= Dy*ur;      % calculate derivative of solution in y
    
  uur = reshape(ur ,mesh_params.ny2,mesh_params.nx)';
  Duur= reshape(Dur,mesh_params.ny2,mesh_params.nx)';
  Dyuur=reshape(Dyur,mesh_params.ny2,mesh_params.nx)';
  
  % Interpolating on the eta mesh
  eta = kz*mesh_params.iz(:) + mesh_params.it(:);  % find skewed rolls and interpolate
  xi  = mesh_params.y(:);
  
  uu_r = zeros(size(eta));
  Duu_r= zeros(size(eta));
  
  uu_r = periodic_interp_2D(uur,mesh_params.x,mesh_params.y2,eta,xi);
  Duu_r= periodic_interp_2D(Duur,mesh_params.x,mesh_params.y2,eta,xi);
  Dyuu_r= periodic_interp_2D(Duur,mesh_params.x,mesh_params.y2,eta,xi);
  
  uu_r = uu_r(:); % reshape for use
  Duu_r= Duu_r(:);
  Dyuu_r=Dyuu_r(:);
end

% interpolate on a periodic grid -vectorised for speed
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
