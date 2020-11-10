%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute 1D Laplacian 4th order finite diffs - r=[0,Lrh]
% Inputs: 2*Lrh - length of domain, Nrh - number of mesh points
% Outputs: r - mesh, L - Laplacian, Dx - 1D differentiation matrix,
% w - integration weights
% Note: returned matrices are sparse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,Dx,D2x,D4x,wx] = Compute_1D_Laplacian_finite_diff_4(nx,Lx)

hx = Lx/(nx-1); x = linspace(0,Lx,nx) - Lx/2;

ex = ones(nx,1); 
Dx = spdiags([ex -8*ex 0*ex 8*ex -ex],-2:2, nx, nx);
Dx(1,:)   = 0; Dx(2,2)   = 1;
Dx(nx,:) = 0; Dx(nx-1,nx-1) = -1;
Dx = Dx/(12*hx);

D2x = spdiags([-ex 16*ex -30.*ex 16*ex -ex], -2:2, nx, nx);
D2x(1,2)= 32; D2x(1,3)= -2;
D2x(2,1)= 16; D2x(2,2)= -31; 
D2x(nx,:) = D2x(1,nx:-1:1);
D2x(nx-1,:) = D2x(2,nx:-1:1);
D2x = D2x/(12*hx^2);

D4x = spdiags([-ex 12*ex -39*ex 56*ex -39*ex +12*ex -ex],-3:3, nx, nx);
D4x(1,2) = -78; D4x(1,3) = 24; D4x(1,4) = -2;
D4x(2,2) =  68; D4x(2,3) = -40;
D4x(3,2) = -40;
D4x(nx,:) = D4x(1,nx:-1:1);
D4x(nx-1,:) = D4x(2,nx:-1:1);
D4x(nx-2,:) = D4x(3,nx:-1:1);
D4x = D4x/(6*hx^4);

Ix = speye(nx);

if (nargout>4)
    wx = [1, 2*ones(1,nx-2)+2*mod([1:nx-2],2),1]/3;   % Simpson weights for intergration int = w*u un-scaled  must *hx
end