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
  np    = p(7);
  
  % Auxiliary variables
  n = mesh_params.nz*mesh_params.nt*mesh_params.ny;
  w     = u(1:n);
  kzm   = u(n+1); 
  kt    = u(n+2); 
  ky    = u(n+3);
  
  chi_p = 1/2 + 1/2*tanh(m*(mesh_params.zz(:)-d));
  chi_m = chi_p(end:-1:1);
  Dchi_p = mesh_params.DZ * chi_p;
  Dchi_m = mesh_params.DZ * chi_m;

  % find skewed rolls
  [uu_rm,Duu_rm,~] = get_sh_rolls(p,kzm,ky,mesh_params,ur0_left);

   w_out = u(1:n); % + chi_m.*uu_rm; % add the chi_m term if want to see original front
   w_out = reshape(w_out,mesh_params.nt,mesh_params.nz,mesh_params.ny);

for i = 1:mesh_params.nt
 pcolor(squeeze(mesh_params.zz(i,:,:)),squeeze(mesh_params.yy(i,:,:))+2*(i-1)*pi,squeeze(w_out(i,:,:)));
 hold on;
end
  
shading interp; axis equal; drawnow;hold off;


   print -dtiff state.tiff

end
