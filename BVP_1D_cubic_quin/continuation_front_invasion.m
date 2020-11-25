clear all, close all, clc;
addpath('../Utilities/');

global ur0_left; % store far-field roll for quick convergence

% Spatial coordinates: x direction
nz = 400; Lz = 4*pi*10; hz = Lz/(nz-1); z = linspace(0,Lz,nz) - Lz/2;

% Temporal discretisation in t
nt = 20; Lt = pi; ht = 2*pi/nt;  t = ht*(1:nt); t = Lt*(t-pi)/pi;

setup_diff_mats; % set up the differentiation matrices/grid etc.

 % Swift-Hohenberg parameters
nu = 1.25;
mu = -0.1; 

kt = 0.99; % initial guess for kt and kz
kz = 0.99;

d = 40; % cut off function parameter: point at x = -d
m = 1;  % cut off function parameter: steepness parameter 

p(1) = mu;
p(2) = nu;
p(5) = m;
p(6) = d;

% partion functions
chi_p = 1/2 + 1/2*tanh(m*(zz(:)-d));
chi_m = chi_p(end:-1:1);

c1p = 1/2*(1+tanh(zz(:)));
c1m = c1p(end:-1:1);

% get far-field rolls
ur0_left = 2*cos(t)';
[uu1_m,~,ur0_left ] = get_sh_rolls(p,kz,mesh_params,ur0_left); % get far-field rolls

w0 = (c1m.*uu1_m - chi_m.*uu1_m);

% setup template functions
mesh_params.w0 = w0;
mesh_params.w0z= DZ*(w0);
mesh_params.w0t= kt*(DT*w0);
mesh_params.loc_wold = (kt/kz*mesh_params.DZ - kt*mesh_params.DT)*mesh_params.w0;
temp = find(zz(:)<(-d+10)); mesh_params.loc_wold(temp) = 0; % eliminate nonzero elements near the cutoff

u0 = [w0; kz; kt]; % initial guess

figure;pcolor(zz,tt,reshape(w0,nt,nz));title('initial condition'); % plot initial guess
shading interp;axis equal;drawnow;

% converge initial guess
my_rhs = @(u) Front_invasion_SH_phase_condition(u,p,mesh_params,0);
options = optimset('Jacobian','on','Display','iter','MaxIter',1000,'DerivativeCheck','off');
[w_kz_out,fval,exitflag,output,jacobian] = fsolve(my_rhs,u0,options);

figure; hold on; % plot converged solution in (z,t) coordinates
w_out = w_kz_out(1:nz*nt);
surf(zz,tt,reshape(w_out,nt,nz));
surf(zz,tt+2*pi,reshape(w_out,nt,nz));
shading interp; axis equal; 
title('converged solution in (z,t)-space');drawnow;

figure;hold on; % plot converged solution in (x,t) coordinates
pcolor(zz+tt,tt,reshape(w_out,nt,nz));
pcolor(zz+(tt+2*pi),tt+2*pi,reshape(w_out,nt,nz));
shading interp; axis equal; 
title('converged solution in (x,t)-space');drawnow;

% Various function handles for continuation
problemHandle            = @(u,p)  Front_invasion_SH_phase_condition(u,p,mesh_params,0);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures(step,u,p,mesh_params);
computeEigenvaluesHandle = [];%@(u,p) ComputeEigenvalues(u,p,t);
plotSpetcrumHandle       = []; %@(d,p,parentHandle) PlotSpectrum(d,p,parentHandle);
% 
stepperPars.iContPar      = 1;
stepperPars.s0            = -0.1;
stepperPars.sMin          = 1e-8;
stepperPars.sMax          = 0.1;
stepperPars.pMin          = -0.28;
stepperPars.pMax          = 100;
stepperPars.maxSteps      = 1000;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','on',...
                                     'MaxIter',15);
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = 'Run_1';
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.PlotBranchVariableId = [];%2;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 1;

branch = SecantContinuation(problemHandle,w_kz_out,p,stepperPars); % continue
