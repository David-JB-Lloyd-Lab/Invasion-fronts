function F = SolutionMeasures_hex(step,u,p,mesh_params)
  
mu = p(1);
nu  = p(2);

  % Auxiliary variables
  n = mesh_params.nx.*mesh_params.ny;
  A     = u(1:n);

  F = [mesh_params.ww2D*(A.^2), u(end)];

end
