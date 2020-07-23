%% fold continuation
function [F,J] = bif_cont(uu,p,mesh_params,problem_handle)

n = mesh_params.nx*mesh_params.ny;
nn= n+2;
u  = uu(1:n);
vv = uu(1+nn:2*nn);
v  = vv(1:end-2);
pp = p;
% pp(1) = uu(2*nn+1);
pp(7) = uu(2*nn+1);

[Fu,Ju] = problem_handle(uu(1:nn),pp);

F = [Fu; Ju*vv; norm(vv,2)^2-1];

if nargout > 1  
    kx = p(8);
    ky = p(7);
    DX = kx*mesh_params.Dx;
    DY = ky*mesh_params.Dy;
    Juu = [spdiags((2*p(2) - 6*u).*v,0,n,n) DX*v DY*v;
           zeros(1,n) 0 0;
           zeros(1,n) 0 0];
           
    J = [Ju  sparse(nn,nn); 
         Juu Ju];
     
    J(2*nn+1,1+nn:2*nn) = 2*vv';
     
     epsiF = 1e-8;
     p2 = pp; p2(7) = p2(7)+epsiF;
     dF = problem_handle(uu(1:nn),p2);
     J(1:nn,2*nn+1) = (dF - Fu)/epsiF;
     J(1+nn:2*nn-2,2*nn+1) = -2*(mesh_params.Dxx*kx^2 + mesh_params.Dyy*ky^2 + mesh_params.Ixy)*2*ky*mesh_params.Dyy*v;
%      J(1+nn:2*nn,2*nn+1) = vv;
%     dF = bif_cont(uu,p2,mesh_params,problem_handle);
%     J(:,2*nn+1) =  (dF - F)/epsiF;
end
     
     
         