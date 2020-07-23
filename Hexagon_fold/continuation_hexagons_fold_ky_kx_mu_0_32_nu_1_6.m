clear all, close all, clc;

addpath('../Utilities');

% Spatial coordinates: x direction
% Spatial discretisation in y
nx = 20; Lx = pi; hx = 2*pi/nx;  x = hx*(1:nx); x = Lx*(x-pi)/pi;
ny = 20; Ly = pi; hy = 2*pi/ny;  y = hy*(1:ny); y = Ly*(y-pi)/pi;

SetupDiffMats; 

 % Swift-Hohenberg parameters
nu = 1.6;
mu = -0.32;
kx = 0.575;
% ky = sqrt(3)/2; % need to change this 
ky = 0.63;

d = 40;
m = 1;
s = 1;

p(1) = mu;
p(2) = nu;
p(5) = m;
p(6) = d;
p(7) = ky;
p(8) = kx;

% partion functions
u0 = 1.4*(cos(ix/kx) + cos((ix/kx+sqrt(3)*iy/ky)/2) + cos((ix/kx-sqrt(3)*iy/ky)/2));

mesh_params.uT = u0(:);

uu0 = [u0(:);0;0];

% sol = load('Hex_ky_run_1/solution_0000060.mat');
% p = sol.p;
% uu0= sol.u;

my_rhs = @(u) Swift_2D_per_cond_2(u,p,mesh_params);

% options = optimset('Jacobian','on','Display','iter','MaxIter',1000,'DerivativeCheck','off','Algorithm','levenberg-marquardt');
options = optimset('Jacobian','on','Display','iter','MaxIter',20,'DerivativeCheck','off');

[u_out,fval,exitflag,output,jacobian] = fsolve(my_rhs,uu0,options);

w_out = u_out(1:nx*ny);
figure;hold on;
pcolor(ix,iy,reshape(w_out,ny,nx));shading interp;axis equal;title('converged hexagon');drawnow;


[V,D] = eigs(jacobian,1,0.01);
myproblemHandle = @(u,p)  Swift_2D_per_cond_2(u,p,mesh_params);
my_rhs2 = @(u) bif_cont(u,p,mesh_params,myproblemHandle);
uu1 = [u_out; V(:,1); p(7)];
options = optimset('Jacobian','on','Display','iter','MaxIter',100,'DerivativeCheck','off');

[u_out2,fval,exitflag,output,jacobian] = fsolve(my_rhs2,uu1,options);

% Various function handles
problemHandle            = @(u,p)  bif_cont(u,p,mesh_params,myproblemHandle);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures_fold(step,u,p,mesh_params);
computeEigenvaluesHandle = [];%@(u,p) ComputeEigenvalues(u,p,t);
plotSpetcrumHandle       = []; %@(d,p,parentHandle) PlotSpectrum(d,p,parentHandle);
% 
stepperPars.iContPar      = 8;
stepperPars.s0            = 0.01;
stepperPars.sMin          = 1e-8;
stepperPars.sMax          = .1;
stepperPars.pMin          = -2;
stepperPars.pMax          = 2;
stepperPars.maxSteps      = 20000;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','on',...
                                     'MaxIter',50);
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = 'Hex_ky_kx_fold_2par_mu_0_32_nu_1_6_run_4';
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.PlotBranchVariableId = [];%2;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 1;

branch = SecantContinuation(problemHandle,u_out2,p,stepperPars);


