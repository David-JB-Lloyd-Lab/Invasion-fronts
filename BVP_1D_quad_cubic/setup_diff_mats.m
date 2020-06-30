%% setup differentiation matrices

% Differentiation matrix in z: finite differences
ez = ones(nz,1); 
Dz = spdiags([ez -8*ez 0*ez 8*ez -ez],-2:2, nz, nz);
Dz(1,:)   = 0; Dz(2,2)   = 1;
Dz(nz,:) = 0; Dz(nz-1,nz-1) = -1;
Dz = Dz/(12*hz);

D2z = spdiags([-ez 16*ez -30.*ez 16*ez -ez], -2:2, nz, nz);
D2z(1,2)= 32; D2z(1,3)= -2;
D2z(2,1)= 16; D2z(2,2)= -31; 
D2z(nz,:) = D2z(1,nz:-1:1);
D2z(nz-1,:) = D2z(2,nz:-1:1);
D2z = D2z/(12*hz^2);

D4z = spdiags([-ez 12*ez -39*ez 56*ez -39*ez +12*ez -ez],-3:3, nz, nz);
D4z(1,2) = -78; D4z(1,3) = 24; D4z(1,4) = -2;
D4z(2,2) =  68; D4z(2,3) = -40;
D4z(3,2) = -40;
D4z(nz,:) = D4z(1,nz:-1:1);
D4z(nz-1,:) = D4z(2,nz:-1:1);
D4z(nz-2,:) = D4z(3,nz:-1:1);
D4z = D4z/(6*hz^4);

Iz = speye(nz);

wz = [1, 2*ones(1,nz-2)+2*mod([1:nz-2],2),1]/3;   % Simpson weights for intergration int = w*u un-scaled  must *hx

% Fourier differentiation matrix first order for t between -pi and pi
column = [0 .5*(-1).^(1:nt-1).*cot((1:nt-1)*ht/2)]';
Dt = toeplitz(column,column([1 nt:-1:2]));

% Fourier differentiation matrix for t between -pi and pi
D2t = toeplitz([-pi^2/(3*ht^2)-1/6 .5*(-1).^(2:nt)./sin(ht*(1:nt-1)/2).^2]);
D4t = D2t^2;

wt = 2*pi*ones(nt,1)/nt; % integration weights for trapzoid rule - mean

wzt = kron(wz,wt'); wzt=wzt(:)';

% Rewrite differentiation matrix for y between 0 and pi
It = speye(nt); 

% Linear 2D differentiation matrices
DT    = kron(Iz,Dt);
DZ    = kron(Dz  ,It);
D2Z   = kron(D2z  ,It);
D2T   = kron(Iz  ,D2t);
D2Z2T = D2Z*D2T;
D4Z   = kron(D4z,It);
D4T   = kron(Iz,D4t);

[zz,tt] = meshgrid(z,t); % 2D mesh

% place all computational matrices etc into mesh_params structure
mesh_params.nz  = nz;  mesh_params.nt  = nt;
mesh_params.Lz  = Lz;  mesh_params.Lt  = Lt;

mesh_params.zz  = zz;  mesh_params.tt  = tt;

mesh_params.Iz  = Iz;  mesh_params.It  = It;
mesh_params.wz  = wz;  mesh_params.wt  = wt;

mesh_params.Dz  = Dz; mesh_params.D2z = D2z; mesh_params.D4z = D4z;

mesh_params.Dt  = Dt; mesh_params.D2t = D2t; mesh_params.D4t = D4t;

mesh_params.D2Z = D2Z; mesh_params.D2T = D2T; mesh_params.DZ = DZ; mesh_params.DT = DT;
mesh_params.D2Z2T = D2Z2T;
mesh_params.D4Z = D4Z; mesh_params.D4T = D4T;
mesh_params.wzt = wzt;
mesh_params.iend=find(abs(zz)==Lz/2);
