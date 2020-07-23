% Continue hexagon eckhaus stability boundary

clear all, close all, clc;

addpath('../Utilities/');

% Spatial coordinates: x and y direction
nx = 16; Lx = pi; hx = 2*pi/nx;  x = hx*(1:nx); x = Lx*(x-pi)/pi;
ny = 16; Ly = pi; hy = 2*pi/ny;  y = hy*(1:ny); y = Ly*(y-pi)/pi;

SetupDiffMats; 

 % Swift-Hohenberg parameters
nu = 1.6;
mu = -0.15;
kx = 0.5;
ky = 0.866; % need to change this

p(1) = mu;
p(2) = nu;
p(3) = ky;
p(4) = kx;

% partion functions
u0 = 3*(cos(ix/kx) + cos((ix/kx+sqrt(3)*iy/ky)/2) + cos((ix/kx-sqrt(3)*iy/ky)/2));

mesh_params.uT = u0(:);

uu0 = [u0(:);0;0];

my_rhs = @(u) Swift_2D_per_cond_2(u,p,mesh_params);

% options = optimset('Jacobian','on','Display','iter','MaxIter',1000,'DerivativeCheck','off','Algorithm','levenberg-marquardt');
options = optimset('Jacobian','on','Display','iter','MaxIter',20,'DerivativeCheck','off');

[u_out,fval,exitflag,output,jacobian] = fsolve(my_rhs,uu0,options);

kx = 0.55; % initial guess for where the eckhaus stability boundary is
ky = 0.65;

p(3) = ky;
p(4) = kx;

my_rhs = @(u) Swift_2D_per_cond_2(u,p,mesh_params);
[u1,fval,exitflag,output,jacobian] = fsolve(my_rhs,u_out,options);

u_out = u1;
w_out = u_out(1:nx*ny);
figure;hold on;
pcolor(ix,iy,reshape(w_out,ny,nx));shading interp;axis equal;drawnow;

% find the eckhaus stability boundary
u1 = [u_out(1:nx*ny);  mesh_params.Dx*u_out(1:nx*ny); 0*ones(nx*ny,1); 0*ones(nx*ny,1); 0; 0; 0; 0; 0];
% 
my_eig = @(u) SH_2D_eig(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out2,fval,exitflag,output,jacobian] = fsolve(my_eig,u1,options);
% 
%
u2 = u_out2;
u2(end) = p(4);
p(1) = 0;
my_eig = @(u) SH_2D_eig_2(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',10,'Algorithm','levenberg-marquardt');
[u_out3,fval,exitflag,output,jacobian] = fsolve(my_eig,u2,options); % find eckhaus boundary


%% Various function handles
problemHandle            = @(u,p)  SH_2D_eig_2(u,p,mesh_params);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution_hex(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures_hex(step,u,p,mesh_params);
computeEigenvaluesHandle = [];%@(u,p) ComputeEigenvalues(u,p,t);
plotSpetcrumHandle       = []; %@(d,p,parentHandle) PlotSpectrum(d,p,parentHandle);
% 
stepperPars.iContPar      = 1;
stepperPars.s0            = -0.1;
stepperPars.sMin          = 1e-8;
stepperPars.sMax          = 0.1;
stepperPars.pMin          = -1.0;
stepperPars.pMax          = 0;
stepperPars.maxSteps      = 20000;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','iter',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','off',...
                                     'MaxIter',15, ...
                                     'Algorithm','levenberg-marquardt');
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = 'Hex_Eckhaus_nu_1_6_ky_0_6_run_1';
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 2;

branch = SecantContinuation(problemHandle,u_out3,p,stepperPars); % continue!

