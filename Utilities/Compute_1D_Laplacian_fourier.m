%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute 1D fourier Laplacian - r=[-Lrh,Lrh]
% Inputs: 2*Lrh - length of domain, Nrh - number of mesh points
% Outputs: r - mesh, L - Laplacian, Dx - 1D differentiation matrix,
% w - integration weights
% Note: returned matrices are sparse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r,L,Dx,w] = Compute_1D_Laplacian_fourier(Nrh,Lrh)

[t,Dt]  = fourdif(Nrh,1); 
[t,Dtt] = fourdif(Nrh,2);



r = Lrh*(t - pi)/pi;

%% 1D
Dx = (pi/Lrh)*Dt;
L  = (pi/Lrh)^2*Dtt;

if (nargout>2)
   w = 2*Lrh*ones(1,Nrh)/Nrh;   % trapezoidal weights for intergration int = w*u
end
