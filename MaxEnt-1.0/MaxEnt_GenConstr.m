
% -----------------------------------------------------------------
%  MaxEnt_GenConstr.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% ----------------------------------------------------------------- 
%  This function computes the MaxEnt distribution for a given set
%  of N independent statistical properties. This distribution is 
%  parametrized by N hyper-parameters defined in terms of the 
%  Lagrange multipliers associated to the underlying optimization
%  convex problem. It employs a robust algorithm for calculation 
%  of the Lagrange multipliers, presented in Soize 2017, which is
%  based on Newton iteration combined with SDV decomposition to
%  solve the underlying linear system (generally ill-conditioned).
%  The algorithm implementation is based on the MaxEnt code presented
%  by Mohammad-Djafari 1992.
%
%  References:
%  C. Soize,
%  Uncertainty Quantification: An Accelerated Course with 
%  Advanced Applications in Computational Engineering.
%  Springer, 2017, pp 221-233
%  
%  A. Mohammad-Djafari,
%  A Matlab Program to Calculate the Maximum Entropy Distributions.
%  In: Smith C.R., Erickson G.J., Neudorfer P.O. (eds)
%  Maximum Entropy and Bayesian Methods, pp 90-92
%  Springer, 1992
%
%  input:
%  xmin    - random variable support left extreme
%  xmax    - random variable support right extreme
%  Nx      - number of points for discretization of X support
%  b       - (N x 1) statistical properties values
%  gfunc   - (N x 1) statistical properties functions
%  lambda0 - (N x 1) Lagrange multipliers initial guess (optional)
%  alpha   - under relaxation factor                    (optional)
%  tol     - numerical tolerance                        (optional)
%  maxiter - maximum number of iterations               (optional)
%
%  output:
%  lambda  - (N  x 1) Lagrange multipliers vector
%  Xpdf    - (Nx x 1) MaxEnt PDF
%  Xsupp   - (Nx x 1) MaxEnt PDF support
%  Xcdf    - (Nx x 1) MaxEnt CDF
%  Xcdfinv - (Nx x 1) MaxEnt quantile function
%  Xprob   - (Nx x 1) MaxEnt quantile function support
%  Entropy - MaxEnt PDF entropy
%  Area    - MaxEnt PDF area
% ----------------------------------------------------------------- 

% -----------------------------------------------------------------
function [lambda,Xpdf,Xsupp,Xcdf,Xcdfinv,Xprob,Entropy,Area] = ...
          MaxEnt_GenConstr(xmin,xmax,Nx,b,gfunc,lambda0,alpha,tol,maxiter)

    % check number of arguments
    if nargin < 5
        error('Too few inputs.')
    elseif nargin > 9
        error('Too many inputs.')
    end
    
    % check arguments
    if xmin >= xmax
        error('xmin must be less than supp_max.')
    end
    
    if Nx < 2
        error('Nx must be an integer greater than one.')
    end
    
    % number of constraints
    N = length(b);
    
    if nargin == 5
        % prealocate memory for lambda0
        lambda0 = zeros(N,1);
        % initial guess for lambda
        lambda0(1) = log(xmax-xmin);
        % under relaxation factor
        alpha = 1.0;
        % tolerance
        tol = 1.0e-6;
        % maximum of iteration
        maxiter = 20;
    elseif nargin == 6
        % under relaxation factor
        alpha = 1.0;
        % tolerance
        tol = 1.0e-6;
        % maximum of iteration
        maxiter = 20;
    elseif nargin == 7
        % tolerance
        tol = 1.0e-6;
        % maximum of iteration
        maxiter = 20;
    elseif nargin == 8
        % maximum of iteration
        maxiter = 20;
    end
    
    % check arguments
    if length(b) ~= length(lambda0)
        error('b and lambda0 vectors must be same length')
    end

    if alpha > 1.0
        error('alpha must be less than one')
    end

    if tol < 0.0
        error('tol must be positive')
    end

    if maxiter <= 1
        error('maxiter must be an integer greather than one')
    end

    % random variable support discretization
    Xsupp = linspace(xmin,xmax,Nx)';

    % quantile function support
    Xprob = linspace(0,1,Nx)';
    
    % preallocate memory for MaxEnt CDF
    Xcdf = zeros(Nx,1);
    
    % statistical properties vector
    G = zeros(N,1);

    % statistical properties functions
	g = gfunc(Xsupp);
    
    % initialize lambda
    lambda = lambda0;
    
    % initialize iteration counter
    iter = 0;
    
    % initialize error estimator
    Error = inf;

    % iteration loop of Newton method
    while Error > tol && iter < maxiter
        
        % update iteration counter
        iter = iter + 1;
        
        % MaxEnt PDF
        Xpdf = exp(-g(:,1:N)*lambda0);
        
        % statistical properties vector
        for n = 1:N
            G(n) = trapz(Xsupp,g(:,n).*Xpdf);
        end
        
        % preallocate memory for statistical properties gradient
        gradG = zeros(N,N);
        
        % statistical properties gradient lower triangle
        gradG(:,1) = -G;
        for m = 2:N
            for n = 2:m
                gradG(m,n) = -trapz(Xsupp,g(:,m).*g(:,n).*Xpdf);
            end
        end
        
        % statistical properties gradient upper triangle
        gradG = gradG + gradG' - diag(diag(gradG));
        
        % right hand side of Newton iteration linear system
        G = b - G;
        
        % SVD factorization of Newton iteration matrix
        [U,sigma,V] = svd(gradG);
        
        % compute lambda correction via SVD
        dlambda = V*((U'*G)./diag(sigma));
        
        % update lambda
        lambda = lambda0 + alpha*dlambda;
        
        % update initial guess
        lambda0 = lambda;
        
        % update error estimator
        Error = norm(G)/norm(b);

    end

    % check for convergence
    if Error > tol
        disp(' ')
        disp('Error estimation in the last iteration:')
        disp(Error)
        disp(' ')
        disp('Approximation obtained in the last iteration:')
        disp(lambda)
        disp(' ')
        error(['Newton method did not converge. ',...
               'Perphaps you initial guess is very far from a solution!'])
    end
    
    % MaxEnt PDF
    Xpdf = exp(-g(:,1:N)*lambda);
    
    % MaxEnt CDF
    for n=2:Nx
        Xcdf(n,1) = trapz(Xsupp(1:n,1),Xpdf(1:n,1));
    end
    
    % quantile function
    Xcdfinv = interp1(Xcdf,Xsupp,Xprob,'linear','extrap')';
    
    % MaxEnt PDF entropy
    Entropy = lambda'*(b-G);
    
    % MaxEnt PDF area
    Area = trapz(Xsupp,Xpdf);
end
% -----------------------------------------------------------------