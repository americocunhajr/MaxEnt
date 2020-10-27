
% -----------------------------------------------------------------
%  MaxEnt_BetaGen.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% ----------------------------------------------------------------- 
%  This functions computes the MaxEnt distribution for the case
%  where the support, mean, the geometric means at the support
%  bounds are the known statistical information i.e., a generalized
%  beta distribution characterized in terms of its support, mean and 
%  standard deviation.
%
%  input:
%  xmin  - support lower bound
%  xmax  - support upper bound
%  mu1   - mean value
%  sigma - standard deviation
%  Nx    - number of points for support discretization
%
%  output:
%  Xpdf    - (Nx x 1) MaxEnt PDF
%  Xsupp   - (Nx x 1) MaxEnt PDF support
%  Xcdf    - (Nx x 1) MaxEnt CDF
%  Xcdfinv - (Nx x 1) MaxEnt quantile function
%  Xprob   - (Nx x 1) MaxEnt quantile function support
%  Entropy - MaxEnt PDF entropy
%  Area    - MaxEnt PDF area
% ----------------------------------------------------------------- 

% -----------------------------------------------------------------
function [Xpdf,Xsupp,Xcdf,Xcdfinv,Xprob,Entropy,Area] = ...
                      MaxEnt_BetaGen(xmin,xmax,mu1,sigma,Nx)

    % check number of arguments
    if nargin < 5
        error('Too few inputs.')
    elseif nargin > 5
        error('Too many inputs.')
    end
    
    % check for consistency
    if Nx < 2
        error('Nx must be greather than or equal to 2')
    end

    if xmin >= xmax
        error('xmax must be greather than xmin');
    end
    
    if mu1 <= xmin || mu1 >= xmax
        error('mu1 must be in (xmin,xmax) interval');
    end
    
    if sigma < 0.0
        error('sigma must be positive');
    end
    
    % standard beta mean
    mu1 = (mu1 - xmin)/(xmax-xmin);
    
    % standard beta standard deviation
    sigma = sigma/(xmax-xmin);
    
    % shape parameters adition
    nu = mu1*(1-mu1)/(sigma^2)-1;

    % shape parameter 1
    shape1 = mu1*nu;

    % shape parameter 2
    shape2 = (1-mu1)*nu;
                      
    % standard beta support
    Xsupp = linspace(0,1,Nx)';
    
    % standard beta PDF function
	Xpdf = betapdf(Xsupp,shape1,shape2);

	% generalized beta support
    Xsupp = xmin + (xmax-xmin)*Xsupp;
                      
    % generalized beta PDF function
	Xpdf = Xpdf/(xmax-xmin);
    
    % preallocate memory for generalized beta CDF function
    Xcdf = zeros(Nx,1);
    
    % generalized beta CDF function
    for n=2:Nx
        Xcdf(n,1) = trapz(Xsupp(1:n,1),Xpdf(1:n,1));
    end
    
    % quantile function support
    Xprob = linspace(0,1,Nx)';
    
    % quantile function
    Xcdfinv = interp1(Xcdf,Xsupp,Xprob,'linear','extrap')';
    
    % Entropy
    Entropy = - trapz(Xsupp,Xpdf.*log(abs(eps+Xpdf)));
    
    % PDF area
    Area = trapz(Xsupp,Xpdf);
end
% -----------------------------------------------------------------
