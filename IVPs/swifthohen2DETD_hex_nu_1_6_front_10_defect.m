% swifthohen2DETD_hex_nu_1_6_front_10_defect.m 
% evolve hexagon <10> front that develops a defect in pattern 
mu = -0.05;
nu = 1.6;
Nx = 2^11; Lx = 100*pi; Ny = 2^6; Ly = 2*pi/0.72;

dt=0.1;                % mesh and domain size and time step
dx=2*pi/Nx; x=dx*(1:Nx )';x=Lx*(x-pi)/pi;
dy=2*pi/Ny; y=dy*(1:Ny )';y=Ly*(y-pi)/pi;
[xx,yy]=meshgrid(x,y);

% Initial data - hexagon front
d = 4*pi;
u = -(-tanh(xx+d)+tanh(xx-d))/2.*((cos(xx) + cos((xx+sqrt(3)*yy)/2) + cos((xx-sqrt(3)*yy)/2)))/2;

% compute eigenvalues of linear operator
k1 = [0:Nx/2 -Nx/2+1:-1]*(pi/Lx);           % wave numbers
k2 = [0:Ny/2 -Ny/2+1:-1]*(pi/Ly);
[kkx,kky]=meshgrid(k1,k2);
LU = - (-(kkx.^2 + kky.^2)+1).^2 + mu;  % Linear Operator
EXP = exp(LU*dt);                       % Exact linear bit
ETD = (exp(LU*dt)-1)./LU;               % ETD1 coeffs of linear operator

uT = fft2(u);
figure;                                 % Plot initial data
pcolor(xx,yy,u);title('t=0');axis equal; axis tight;shading interp;drawnow

 for t=0:dt:300                        % evolve via ETD  method
     f = nu*u.^2 - u.^3;                 % nonlinear rhs
     
     % exponential timestep in fourier space
     fT = fft2(f); uT = uT.*EXP + fT.*ETD;      
     u   = real(ifft2(uT));                         

     if mod(t,1) == 0                 % display at intervals
         pcolor(xx,yy,u);%view([0 90]);
         title(['t=' num2str(t) ' max u=' num2str(max(max(u)))]);
         shading interp; axis equal;axis tight; drawnow;
     end
 end
 