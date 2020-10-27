
% -----------------------------------------------------------------
%  MaxEnt_Uniform.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% ----------------------------------------------------------------- 
%  This functions computes the MaxEnt distribution for the case
%  where the support is the only known statistical information,
%  i.e., an uniform distribution given the support.
%
%  input:
%  xmin  - support lower bound
%  xmax  - support upper bound
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
                                    MaxEnt_Uniform(xmin,xmax,Nx)

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
    
    if xmin >= xmax
        error('xmin must be less than xmax');
    end
                              
    % PDF support
    Xsupp = linspace(xmin,xmax,Nx)';
    
    % PDF function
	Xpdf = unifpdf(Xsupp,xmin,xmax);
    
    % CDF function
    Xcdf = unifcdf(Xsupp,xmin,xmax);
    
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
