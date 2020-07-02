% Main continuation script for continuation of traveling <10> hexagon front
% Solves the Swift-Hohenberg equation with a non-variational term
% u_t = -(1+Delta)^2u + mu*u - nu*u^2 - u^3 + alpha*|nabla*u|^2 term

clear all, close all, clc;
addpath('../Utilities');

global ur0_left; % make global previous far-field hexagon cell. 

% Spatial coordinates: x direction
nz = 400; Lz = 5*pi*10; Lend = 10*pi; z = linspace(-Lz/2,Lend,nz); hz = z(2)-z(1);

% Spatial discretisation in y
ny = 16; Ly = pi; hy = 2*pi/ny;  y = hy*(1:ny); y = Ly*(y-pi)/pi;

% spatial disceretisation for a single hexagon
nx = 16; Lx = pi; hx = 2*pi/nx;  x = hx*(1:nx); x = Lx*(x-pi)/pi;
ny2= 16; Ly2= pi; hy2= 2*pi/ny2; y2= hy2*(1:ny2);y2 = Ly2*(y2-pi)/pi;

% Temporal discretisation in t
nt = 16; Lt = pi; ht = 2*pi/nt;  t = ht*(1:nt); t = Lt*(t-pi)/pi;

SetupDiffMats; % build differentiation matrices, mesh etc. 

 % Swift-Hohenberg parameters
mu = -0.1;
nu = 1.6;
kt = 2;     % initial guess wavenumber in time
ky = 0.75;  % wavenumber in y
kz = 0.5;   % initial guess wavenumber in z

d = 40;     % cutoff function distance at z = -d
m = 1;      % scaling of shape of cutoff function
s = 1;

alpha = 0; % non-var alpha*|nabla*u|^2 term

p(1) = mu;      % set up p vector for continuation
p(2) = nu;
p(5) = m;
p(6) = d;
p(7) = ky;
p(8) = alpha; 

% partion functions
chi_m = 1/2 + 1/2*tanh(m*(-zz(:)-d));
c1m = 1/2*(1+tanh(-zz(:)));

% setup initial condition for far-field single hexagon cell
ur0_left = 1.4/3*(cos(ix/kz) + cos((ix/kz+sqrt(3)*iy/ky)/2) + cos((ix/kz-sqrt(3)*iy/ky)/2));
ur0_left = ur0_left(:); mesh_params.uT = ur0_left;
[uu1_m,~,ur0_left ] = get_sh_hexs_nonvar(p,kz,ky,mesh_params,ur0_left); % get/converge far-field hexagon

w0 = (c1m.*uu1_m - chi_m.*uu1_m); % for initial traveling front guess

w_out = reshape(w0(1:nz*nt*ny),nt,nz,ny); % plot initial guess
figure;hold on;
for i = 1:nt
pcolor(squeeze(kz*zz(i,:,:)+tt(i,:,:)+pi),squeeze(yy(i,:,:))+2*(i-1)*pi,squeeze(w_out(i,:,:)));
axis equal;shading interp;drawnow
end

mesh_params.w0 = w0;        % set-up key template functions
mesh_params.w0z= DZ*(w0);
mesh_params.w0t= kt*(DT*w0);
mesh_params.loc= (DZ*kt/kz - kt*DT)*w0;
temp = find(zz(:)<=-d+10);mesh_params.loc(temp) = 0;

u0 = [w0; kz; kt];          % initial guess to be continued

% converge initial front guess
my_rhs = @(u) Front_invasion_SH_phase_condition_nonvar(u,p,mesh_params,0); 
options.nonlinTol = 1e-4;
options.nonlinMaxIter = 20;
options.display = 1;
options.border = nz*nt*ny;

[w_kz_out,resHist,flag] = Newton_Amijo(my_rhs,u0,options); 

% plot converged solution
w_out = reshape(w_kz_out(1:nz*nt*ny),nt,nz,ny);
figure;hold on;
for i = 1:nt
pcolor(squeeze(w_kz_out(end-1)*zz(i,:,:)+tt(i,:,:)+pi),squeeze(yy(i,:,:))+2*(i-1)*pi,squeeze(w_out(i,:,:)));
axis equal;shading interp;
end
drawnow;

% continue front in mu
problemHandle            = @(u,p)  Front_invasion_SH_phase_condition_nonvar(u,p,mesh_params,0);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures(step,u,p,mesh_params);
computeEigenvaluesHandle = [];%@(u,p) ComputeEigenvalues(u,p,t);
plotSpetcrumHandle       = []; %@(d,p,parentHandle) PlotSpectrum(d,p,parentHandle);
stepperPars.iContPar      = 1; % continue in p(iContPar)
stepperPars.s0            = 1; % initial s0 guess +/- to change direction
stepperPars.sMin          = 1e-8; % min arclength step
stepperPars.sMax          = 2;    % max arclength step
stepperPars.pMin          = -1;   % stop if p(iContPar)<pMin
stepperPars.pMax          = 0; % stop if p(iContPar)>pMax
stepperPars.maxSteps      = 1000; % max number of continuation steps
stepperPars.nPrint        = 1;    % Print solution at every step
stepperPars.nSaveSol      = 1;    % save solution at every step
stepperPars.finDiffEps    = 1e-7; % finite-difference step for pseudo-arclength continuation
stepperPars.fsolveOptions = options; % my_fsolve border options
stepperPars.optNonlinIter = 10; % not used
stepperPars.dataFolder    = 'Hex_10_non_var_fixed_ky_0_75_Nt_16_Ny_16_Nz_400_run_1';
stepperPars.PlotSolution  = plotSolutionHandle; % plot solution function
stepperPars.BranchVariables = branchVariablesHandle; % branch solution function
stepperPars.PlotBranchVariableId = [];%2;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 1; % plot the solutionMeasures ith variable

% continue!
branch = SecantContinuation_bordered_Amijo(problemHandle,w_kz_out,p,stepperPars);
