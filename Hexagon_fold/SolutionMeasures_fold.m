function F = SolutionMeasures(step,u,p,mesh_params)
  
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
  n = mesh_params.nx*mesh_params.ny;
  w     = u(1:n);
  cy    = u(n+1);
  cx    = u(n+2);

   w_out = u(1:n); 
   
F = [u(end) cx cy];
  
  

end
