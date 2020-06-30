function plotHandle = PlotSolution(u,p,parentHandle,mesh_params)

global ur0_left;

   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[2/4*scrsz(3) scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
   else
     plotHandle = parentHandle;
   end 

   figure(parentHandle); %hold on;
   
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
 
  %
  chi_p = 1/2 + 1/2*tanh(m*(mesh_params.zz(:)-d));
  chi_m = chi_p(end:-1:1);
  Dchi_p = mesh_params.DZ * chi_p;
  Dchi_m = mesh_params.DZ * chi_m;

  % find skewed rolls
  [uu_rm,Duu_rm,~] = get_sh_rolls(p,kzm,mesh_params,ur0_left);
   
   nz = mesh_params.nz;
   nt = mesh_params.nt;
   w_out = u(1:nz*nt); + chi_m.*uu_rm;
   w_out = reshape(w_out,nt,nz);
 
    pcolor([mesh_params.zz; mesh_params.zz(1,:)],...
           [mesh_params.tt/kt; (mesh_params.tt(1,1)+2*pi)*ones(size(mesh_params.zz(1,:)))/kt],...
 	  [w_out; w_out(1,:)]);
  
shading interp; axis equal; drawnow;hold off;


   print -dtiff state.tiff

end
