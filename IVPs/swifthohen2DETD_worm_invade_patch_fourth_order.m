% swifthohen2DETD_worm_invad e_patch_fourth_order.m
% evolve a strip patch in the cubic-quintic Swift-Hohenberg equation
% code is modified from:
% Kassam, A.-K. (2003). Solving reaction-diffusion equations 10 times faster.
% https://ora.ox.ac.uk/objects/uuid:77b28856-b1fb-418a-a252-65c80cb62d3d

mu = -0.2; % equation parameters
nu = 1.25;

dt=0.01;                % time step
N = 2^10;               % number of Fourier modes in x and y
L = 35*pi;              % (x,y) in [-L,L]^2
dx=2*pi/N;x=dx*(1:N)';x=L*(x-pi)/pi;y=x'; % mesh
[xx,yy]=meshgrid(x,y);

% Initial data - stripe patch
worm = cos(xx); + 0.1*cos(yy);
d1 = 8*pi;
d2 = 4*pi;
u = (-tanh(xx-d1) + tanh(xx+d1)).*(-tanh(yy-d2) + tanh(yy+d2)).*worm/4*1.2;

% compute eigenvalues of linear operator
k = [0:N/2 -N/2+1:-1]*(pi/L);           % wave numbers
[kkx,kky]=meshgrid(k,k);
LU = - (-(kkx.^2 + kky.^2)+1).^2 + mu;   % Linear Operator

Fr=logical(zeros(N,1)); %High frequencies for de-aliasing 
Fr([N/2+1-round(N/6) : N/2+round(N/6)])=1; 
[alxi,aleta]=ndgrid(Fr,Fr);
ind=alxi | aleta; 

%=============== PRECOMPUTING ETDRK4 COEFFS ===================== 
E=exp(dt*LU); E2=exp(dt*LU/2);
M=16; % no. of points for complex mean
r=exp(1i*pi*((1:M)-0.5)/M); % roots of unity
LU=LU(:); LR=dt*LU(:,ones(M,1))+r(ones(N^2,1),:); 
Q=dt*real(mean( (exp(LR/2)-1)./LR ,2)); 
f1=dt*real(mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2)); 
f2=dt*real(mean( (4+2*LR+exp(LR).*(-4+2*LR))./LR.^3 ,2)); 
f3=dt*real(mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2)); 
f1=reshape(f1,N,N); f2=reshape(f2,N,N); f3=reshape(f3,N,N); 
LU=reshape(LU,N,N); Q=reshape(Q,N,N); clear LR 

% Locate edges of the patch at midpoints in x and y (u=0.5)
F1 = griddedInterpolant(x,abs(u(N/2,:)));
F2 = griddedInterpolant(y,abs(u(:,N/2)));

iix=find(abs(diff(sign(F1(x)-0.5)))==2);
iiy=find(abs(diff(sign(F2(y)-0.5)))==2);
yfront = fzero(@(x) F1(x)-0.5,x(iix(1)));
xfront = fzero(@(x) F2(x)-0.5,y(iiy(1)));

% FFT intial condition and set up nonlinear function
v = fft2(u);
g=inline('nu*u.^3 - u.^5','u','nu'); %nonlinear function

h = figure(1);                                 % Plot initial data
set(h,'color','w');
set(h,'nextplot','replacechildren');
subplot(1, 2, 1);
pcolor(xx,yy,u);title('t=0');shading interp;axis equal;axis tight;
hold on;
plot(-yfront,0,'or','LineWidth',2);
plot(0,-xfront,'om','LineWidth',2);
hold off;
subplot(1, 2, 2); 
hold on;
plot(0,-xfront,'om','LineWidth',2);
plot(0,-yfront,'or','LineWidth',2);
hold off;axis square;axis([0 50 0 90]); 
drawnow;

Tend = 100; % time step to T = Tend
tt  = 0;
 for t=dt:dt:100                      % evolve via ETD  method
     tt = [t; tt];
     Nv=fftn( g(real(ifftn(v)),nu) ); 
     a=E2.*v + Q.*Nv;
     Na=fftn( g(real(ifftn(a)),nu) ); 
     b=E2.*v + Q.*Na;
     Nb=fftn( g(real(ifftn(b)),nu) );
     c=E2.*a + Q.*(2*Nb-Nv);
     Nc=fftn( g(real(ifftn(c)),nu) );
     v=E.*v + Nv.*f1 + (Na+Nb).*f2 + Nc.*f3; %update 
     v(ind) = 0; % High frequency removal --- de-aliasing
  
     u   = real(ifft2(v));     
     
     % track location of midpoint front interfaces
     F1 = griddedInterpolant(x,abs(u(N/2,:)));
     F2 = griddedInterpolant(y,abs(u(:,N/2)));

     iiy=find(abs(diff(sign(F1(x)-0.5)))==2);
     iix=find(abs(diff(sign(F2(y)-0.5)))==2);
     yfront = [yfront; fzero(@(x) F1(x)-0.5,x(iiy(1)))];
     xfront = [xfront; fzero(@(x) F2(x)-0.5,y(iix(1)))]; 
            
     if mod(t,10*dt) == 0                 % display at intervals
         subplot(1, 2, 1);
         pcolor(xx,yy,u);title('t=0');shading interp;axis equal;axis tight;
         title(['t=' num2str(t) ' max u=' num2str(max(max(u)))]);
         hold on;
         plot(-yfront(end),0,'or','LineWidth',2);
         plot(0,-xfront(end),'om','LineWidth',2);
         hold off;
         subplot(1, 2, 2); 
         
         plot(tt,-xfront,'b','LineWidth',2);hold on;
         plot(tt,-yfront,'b','LineWidth',2);
         plot(t,-xfront(end),'om','LineWidth',2);
         plot(t,-yfront(end),'or','LineWidth',2);
         hold off;axis square;axis([0 100 0 90]); 
         drawnow;
         
     end
 end
TT = 0:dt:(Tend);
figure;plot(TT,-xfront,'b',TT,-yfront,'r');
axis square; xlabel('time'); ylabel('patch interface locations'); title('Patch midpoint interface locations');
legend('xfont','yfront');

figure;
plot(k,abs(fft(u(N/2,:))));title('FFT of midpoint in y');xlabel('k');ylabel('abs(fft(u(:,N/2)))');
