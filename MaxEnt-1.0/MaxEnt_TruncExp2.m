
% -----------------------------------------------------------------
%  MaxEnt_TruncExp2.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% ----------------------------------------------------------------- 
%  This functions computes the MaxEnt distribution for the case
%  where the support and mean are the known statistical information
%  i.e., a truncted exponential characterized in terms of two 
%  parameters (Lagrange multipliers).
%
%  input:
%  xmin   - support lower bound
%  xmax   - support upper bound
%  lambda - (2 x 1) Lagrange multipliers vector
%  Nx     - number of points for support discretization
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
                                    MaxEnt_TruncExp2(xmin,xmax,lambda,Nx)

    % check number of arguments
    if nargin < 4
        error('Too few inputs.')
    elseif nargin > 4
        error('Too many inputs.')
    end
    
    % check for consistency
    if Nx < 2
        error('Nx must be greather than or equal to 2')
    end
    
    if xmin >= xmax
        error('xmin must be less than xmax');
    end
    
    % ensure lambda is a column vector
    lambda = lambda(:);
    
    % check if lambda is a 2-dimensional vector
    if length(lambda) ~= 2
        error('lambda must be a (2 x 1) array')
    end
        
	% Lagrange multipliers
	L0 = lambda(1);
	L1 = lambda(2);
    
    % check consistency
    if L1 <= eps
        error('lambda1 must be a nonzero parameter')
    end
    
    % PDF support
    Xsupp = linspace(xmin,xmax,Nx)';
        
	% PDF function
	Xpdf = exp(-L0-L1*Xsupp);
    
	% CDF function
	Xcdf = exp(-L0)*(exp(-L1*xmin)-exp(-L1*Xsupp))/L1;
    
    % quantile function support
    Xprob = linspace(0,1,Nx)';
    
	% quantile function
	%Xcdfinv = -log(abs(exp(- L1*xmin) - L1*exp(L0)*Xsupp))/L1;
    Xcdfinv = interp1(Xcdf,Xsupp,Xprob,'linear','extrap')';
    
    % Entropy
    Entropy = - trapz(Xsupp,Xpdf.*log(abs(eps+Xpdf)));
    
    % PDF area
    Area = trapz(Xsupp,Xpdf);

return
% -----------------------------------------------------------------