% Main continuation script for continuation of invading stripe fronts in the
% cubic-qunitic Swift-Hoheberg equation
clear all, close all, clc;

global ur0_left;    % make global previous far-field stripe
addpath('../Utilities/');

% Spatial coordinates: x direction
nz = 400; Lz = 4*pi*10; hz = Lz/(nz-1); z = linspace(0,Lz,nz) - Lz/2;

% Spatial discretisation in y
nx = 16; Lx = pi; hx = 2*pi/nx;  x = hx*(1:nx); x = Lx*(x-pi)/pi;
ny = 16; Ly = pi; hy = 2*pi/ny;  y = hy*(1:ny); y = Ly*(y-pi)/pi;

% Temporal discretisation in t
nt = 16; Lt = pi; ht = 2*pi/nt;  t = ht*(1:nt); t = Lt*(t-pi)/pi;

SetupDiffMats; % build differentiation matrices, mesh etc.

 % Swift-Hohenberg parameters
nu = 1.25;
mu = -0.2;
kt = 2;   % initial guess wavenumber in time
ky = 0.5;   % wavenumber in y
kz = 0.995;  % initial guess wavenumber in z

d = 40;     % cutoff function distance at z = -d
m = 1;      % scaling of shape of cutoff function
s = 1;

p(1) = mu;  % set up p vector for continuation
p(2) = nu;
p(5) = m;
p(6) = d;
p(7) = ky;

% partion functions
chi_p = 1/2 + 1/2*tanh(m*(zz(:)-d));
chi_m = chi_p(end:-1:1);

c1p = 1/2*(1+tanh(zz(:)));
c1m = c1p(end:-1:1);

ur0_left = 2*cos(x'/kz);
ur0_left = ur0_left(:);
mesh_params.uT = ur0_left;
[uu1_m,~,ur0_left ] = get_sh_rolls(p,kz,ky,mesh_params,ur0_left); % get far-field rolls

w0 = (c1m.*uu1_m - chi_m.*uu1_m); % for initial traveling front guess
dd = pi/2;
loc= (-tanh(zz(:)-dd)+tanh(zz(:)+dd))/2;
w0 = w0 + loc.*cos(yy(:)/ky).*cos(zz(:)/kz);

mesh_params.w0 = w0;            % set-up key template functions
mesh_params.w0z= DZ*(w0);
mesh_params.w0t= kt*(DT*w0);
mesh_params.loc_wold =  (kt/kz*mesh_params.DZ - kt*mesh_params.DT)*mesh_params.w0;
temp = find(zz(:)<=(-d+10)); mesh_params.loc_wold(temp) = 0;

% u0 = [w0; kz; kt; 0]; % initial guess to be continued - this requires a
% lot of newton iterations so use previously converged solution
sol = load('solution_0000000.mat');
u0  = sol.u;
p   = sol.p;

% converge initial front guess
my_rhs = @(u) Front_invasion_SH_phase_condition(u,p,mesh_params,0);
options.nonlinTol = 1e-4;
options.nonlinMaxIter = 20;
options.display = 1;
options.border = nz*nt*ny;

[w_kz_out,resHist,flag] = Newton_Amijo(my_rhs,u0,options);

% plot converged solution
figure; hold on;
w_out = reshape(w_kz_out(1:nz*nt*ny),nt,nz,ny);
for i =1:nt
   pcolor(squeeze(zz(i,:,:)),squeeze(yy(i,:,:))+2*(i-1)*pi,squeeze(w_out(i,:,:)));
end
shading interp; axis equal; title('Converged'); drawnow;


% continue front in mu
problemHandle            = @(u,p)  Front_invasion_SH_phase_condition(u,p,mesh_params,0);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures(step,u,p,mesh_params);
computeEigenvaluesHandle = [];%@(u,p) ComputeEigenvalues(u,p,t);
plotSpetcrumHandle       = []; %@(d,p,parentHandle) PlotSpectrum(d,p,parentHandle);
% 
stepperPars.iContPar      = 1;      % continue in p(iContPar)
stepperPars.s0            = -0.1;   % initial s0 guess +/- to change direction
stepperPars.sMin          = 1e-8;   % min arclength step
stepperPars.sMax          = 0.1;    % max arclength step
stepperPars.pMin          = -0.6;   % stop if p(iContPar)<pMin
stepperPars.pMax          = 100;    % stop if p(iContPar)>pMax
stepperPars.maxSteps      = 2000;   % max number of continuation steps
stepperPars.nPrint        = 1;      % Print solution at every step
stepperPars.nSaveSol      = 1;      % save solution at every step
stepperPars.finDiffEps    = 1e-7;   % finite-difference step for pseudo-arclength continuation
stepperPars.fsolveOptions = options; % my_fsolve border options
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = 'Run_1_nu_1_25';
stepperPars.PlotSolution  = plotSolutionHandle;         % plot solution function
stepperPars.BranchVariables = branchVariablesHandle;    % branch solution function
stepperPars.PlotBranchVariableId = [];%2;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 1;       % plot the solutionMeasures ith variable

% continue!
branch = SecantContinuation_bordered_Amijo(problemHandle,w_kz_out,p,stepperPars);
