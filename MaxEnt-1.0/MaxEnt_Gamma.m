
% -----------------------------------------------------------------
%  MaxEnt_Gamma.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% ----------------------------------------------------------------- 
%  This functions computes the MaxEnt distribution for the case
%  where the support [0,infty], mean and the geometric mean are
%  the known statistical information i.e., a gamma distribution
%  characterized in terms of its mean, standard deviation and a 
%  finite support for purposes of computation.
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
                            MaxEnt_Gamma(xmin,xmax,mu1,sigma,Nx)

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
        error('xmin must be less than xmax');
    end
    
    if mu1 <= xmin || mu1 >= xmax
        error('mu1 must be in (xmin,xmax) interval');
    end
    
    if sigma < 0.0
        error('sigma must be positive');
    end

    % shape parameter
    shape = (mu1/sigma)^2;

    % scale parameter
    scale = sigma^2/mu1;
                      
    % PDF support
    Xsupp = linspace(xmin,xmax,Nx)';
    
    % PDF function
	Xpdf = gampdf(Xsupp,shape,scale);
    
    % CDF function
    Xcdf = gamcdf(Xsupp,shape,scale);
    
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
