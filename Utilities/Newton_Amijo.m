% Newton_Amijo code from
% Solving Nonlinear Equations with Newton's Method, C. T. Kelley
% Published: 2003
% ISBN: 978-0-89871-546-0
% eISBN: 978-0-89871-889-8
% https://doi.org/10.1137/1.9780898718898
%
% Bordering code changes by David J.B. Lloyd
function [x,resHist,flag,output] = Newton_Amijo(fhandle,x0,options)

  % Rename parameters
  nltol       = options.nonlinTol;
  nlmaxit     = options.nonlinMaxIter;
  display     = options.display;
  border      = options.border; % the amount of sparse structure

  % Initialise iterations
  x = x0;
  [f] = feval(fhandle,x);
  neval = 1;
  res = norm(f,2);
  resHist = res;
  it = 0;

  alpha = 1e-4;
  sigma0 = 0.1;
  sigma1 = 0.5;
  maxarm = 20;
  gamma = 0.9;
  
  res0 = 1;
  % Displaying results
  if display
    tic
    DisplayIteration(it,neval,res);
  end
  
  output.iterations = 0;

  % Main loop
  while ( (res > nltol) && (it < nlmaxit) )

    % Linear solve bordered system   
    [f,J] = feval(fhandle, x);
    
    JJ = J(1:border,1:border);
    A  = J(1:border,1+border:end);
    BT = J(1+border:end,1:border);
    C  = J(1+border:end,1+border:end);
    FF = -f(1:border);
    G  = -f(1+border:end);
    
    d = solve_bordered_sys(JJ,A,BT,C,FF,G);
    
    rat = res/res0;
    res0 = res;
    
    % Update solution

    xold = x;
    lambda = 1; lamm = 1; lamc = lambda; iarm = 0;
    xt = x + lambda*d;
    [ft] = feval(fhandle,xt);
    nft = norm(ft); nf0 = res; ff0 = nf0*nf0; ffc = nft*nft; ffm = ffc;
    while nft >= (1 - alpha*lambda) * nf0
          if iarm == 0
              lambda = sigma1*lambda;
          else
              lambda = parab3p(lamc, lamm, ff0, ffc, ffm);
          end
          
          xt = x + lambda*d;
          lamm = lamc;
          lamc = lambda;
          
          [ft] = feval(fhandle,xt);
          nft = norm(ft);
          ffm = ffc;
          ffc = nft^2;
          iarm = iarm + 1;
          if iarm > maxarm
              warning(['Armijo failure, too many reductions']);break;
          end
    end
    x   = xt;
    res0= ft;

    % Book-keeping
    [f] = feval(fhandle, x); neval = neval + 1;
    res = norm(f,2); resHit = [resHist; res];
    it = it + 1;   output.iterations = it;

    % Displaying results
    if display
      DisplayIteration(it,neval,res);
    end

  end

  % Output flag
  if res > nltol || isnan(res)
    flag = 0; 
  else
    flag = 1; 
  end
  if display
    toc
  end

  function DisplayIteration(i,funceval,residual)
    if i == 0
      message = sprintf(['\n Start Newton-Armijo with bordering Iterations \n' ...
                           '   Iterations      Func-count      f(x)\n']);
      fprintf(message);
    end
    message = sprintf('%9d %16d %14.4e\n', i, funceval, residual);
    fprintf(message);
  end

 function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
        sigma0 = 0.1;
        sigma1 = 0.5;
        
        c2 = lambdam*(ffc - ff0) - lambdac*(ffm - ff0);
        if c2 >= 0
            lambdap = sigma1*lambdac; return
        end
        c1 = lambdac^2*(ffm - ff0) - lambdam^2*(ffc - ff0);
        lambdap = -c1*0.5/c2;
        if lambdap < sigma0*lambdac, lambdap = sigma0*lambdac; end
        if lambdap > sigma1*lambdac, lambdap = sigma1*lambdac; end
 end

%% solve bordered linear system where n >> m 
% i.e. the number of phase conditions is m
% Solves
% [J  A]X = [F]
% [BT C]    [G]
% requires:
% J = n x n matrix, A = n x m matrix, BT= m x n matrix, C = m x m matrix
% F = n vector, G = m vector 
% incomplete LU decomposition of J
% returns:
% [Y; Z] where Y is a vector of length n and Z is a vector of
% length m
function X = solve_bordered_sys(J,A,BT,C,F,G)
 
[n,m] = size(A);
% solve bordered system
Temp1 = [C'; BT'];

[R,W] = house_qr(Temp1); R = R(1:m,:); % elminate the zero part

T = compactWY(W,ones(m,1)); T = -T;  % Q = speye(n+m)+W*T*W'; % This is the same as [Q,R] = qr(Temp1)
W1 = W(1:m,:); W2 = W(m+1:end,:);
U = (A*W1 + J*W2)*T; VT= W2';  

ZY = R'\G; 

F2 = F - [A J]*house_apply(W,[ZY;zeros(n,1)]);

dJ = decomposition(J); % carry out a sparse LU decomposition of J
Z = dJ\U; Y = dJ\F2; % solve for Z and Y using the LU decomposition

ZX= Y - Z*(inv(speye(m)+VT*Z)*VT*Y); % update using Woodbury matrix formula

X = house_apply(W,[ZY;ZX]); % compute [Y;X] = Q*[ZY;ZX] without computing Q!
X = [X(m+1:end); X(1:m)];
end

function [R,U] = house_qr(A)
    % Householder reflections for QR decomposition.
    % [R,U] = house_qr(A) returns
    % R, the upper triangular factor, and
    % U, the reflector generators for use by house_apply.   
    % code from 
    % https://blogs.mathworks.com/cleve/2016/10/03/householder-reflections-and-the-qr-decomposition/
    H = @(u,x) x - u*(u'*x);
    [m,n] = size(A);
    U = zeros(m,n);
    R = A;
    for j = 1:min(m,n)
        u = house_gen(R(j:m,j));
        U(j:m,j) = u;
        R(j:m,j:n) = H(u,R(j:m,j:n));
        R(j+1:m,j) = 0;
    end
end

function u = house_gen(x)
    % u = house_gen(x)
    % Generate Householder reflection.
    % u = house_gen(x) returns u with norm(u) = sqrt(2), and
    % H(u,x) = x - u*(u'*x) = -+ norm(x)*e_1.
    
    % Modify the sign function so that sign(0) = 1.
    sig = @(u) sign(u) + (u==0);
    
    nu = norm(x);
    if nu ~= 0
        u = x/nu;
        u(1) = u(1) + sig(u(1));
        u = u/sqrt(abs(u(1)));
    else
        u = x;
        u(1) = sqrt(2);
    end
end

function Z = house_apply(U,X)
    % Apply Householder reflections.
    % Z = house_apply(U,X), with U from house_qr
    % computes Q*X without actually computing Q.
    % code from
    % https://blogs.mathworks.com/cleve/2016/10/03/householder-reflections-and-the-qr-decomposition/
    H = @(u,x) x - u*(u'*x);
    Z = X;
    [~,n] = size(U);
    for j = n:-1:1
        Z = H(U(:,j),Z);
    end
end

function T = compactWY(V,beta)
% Computes the compact WY representation of a product of k Householder % reflectors.
% Input: V - m-by-k matrix containing the vectors v_j 
%        beta - vector of length k containing the scalar factors beta_j %
%        note these should just be one's
% Output: T - upper triangular factor of the compact WY representation
% code from: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.184.2228&rep=rep1&type=pdf
[m,k] = size(V); T = zeros(k); 
for j = 1:k
    T(1:j-1,j) = -beta(j) * T(1:j-1,1:j-1) * ( V(:,1:j-1)'*V(:,j) );
    T(j,j) = beta(j);
end


end
end
