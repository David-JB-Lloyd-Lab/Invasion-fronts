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
mu = p(1);
nu = p(2);
ky = p(3);
kx = p(4);
  
  % Auxiliary variables
  n = mesh_params.nx*mesh_params.ny;
  w = u(1:n);
 

 
pcolor(mesh_params.ix,mesh_params.iy,reshape(w,mesh_params.ny,mesh_params.nx));
  
shading interp; axis equal; drawnow;hold off;


   print -dtiff state.tiff

end
