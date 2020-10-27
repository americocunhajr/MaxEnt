
% -----------------------------------------------------------------
%  MaxEnt_Beta.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% ----------------------------------------------------------------- 
%  This functions computes the MaxEnt distribution for the case
%  where the support [0,1] and the geometric means at the support
%  bounds are the known statistical information i.e., a standard
%  beta distribution characterized in terms of its mean and 
%  standard deviation.
%
%  input:
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
                      MaxEnt_Beta(mu1,sigma,Nx)

    % check number of arguments
    if nargin < 3
        error('Too few inputs.')
    elseif nargin > 3
        error('Too many inputs.')
    end
    
    % check for consistency
    if Nx < 2
        error('Nx must be greather than or equal to 2')
    end

    % PDF support limits
    xmin = 0.0;
    xmax = 1.0;
    
    if mu1 <= xmin || mu1 >= xmax
        error('mu1 must be in (xmin,xmax) interval');
    end
    
    if sigma < 0.0
        error('sigma must be positive');
    end
    
    % shape parameters adition
    nu = mu1*(1-mu1)/(sigma^2)-1;

    % shape parameter 1
    shape1 = mu1*nu;

    % shape parameter 2
    shape2 = (1-mu1)*nu;
                      
    % PDF support
    Xsupp = linspace(0,1,Nx)';
    
    % PDF function
	Xpdf = betapdf(Xsupp,shape1,shape2);
    
    % CDF function
    Xcdf = betacdf(Xsupp,shape1,shape2);
    
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
