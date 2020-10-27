
% -----------------------------------------------------------------
%  MaxEnt_TruncExp3.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% ----------------------------------------------------------------- 
%  This functions computes the MaxEnt distribution for the case
%  where the support and mean are the known statistical information
%  i.e., a truncted exponential characterized in terms of three
%  parameters (Lagrange multipliers).
%
%  input:
%  xmin   - support lower bound
%  xmax   - support upper bound
%  lambda - (3 x 1) Lagrange multipliers vector
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
                                    MaxEnt_TruncExp3(xmin,xmax,lambda,Nx)

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
    
    % check if lambda is a 3-dimensional vector
    if length(lambda) ~= 3
        error('lambda must be a (3 x 1) array')
    end
        
	% Lagrange multipliers
	L0 = lambda(1);
	L1 = lambda(2);
    L2 = lambda(3);
    
    % check consistency
    if abs(L2) <= eps
        error('lambda2 must be a nonzero parameter')
    end
    
    % PDF support
    Xsupp = linspace(xmin,xmax,Nx)';
    
	% PDF function
	Xpdf = exp(-L0-L1*Xsupp-L2*Xsupp.^2);
    
    % auxiliar variables
    a = exp(-L0 + 0.25*L1^2/L2);
    b = L2;
    c = -0.5*L1/L2;
        
    if b > 0
        
        % auxiliar variables
        A = 0.5*sqrt(pi)*a/sqrt(b);
        xminus = sqrt(b)*(c-xmin);
        xplus  = sqrt(b)*(c-Xsupp);
        
        % CDF function
        Xcdf = A*(erf(xminus)-erf(xplus));
        
        % quantile function
        %Xcdfinv = c - erfinv(-Xsupp/A + erf(xminus))/sqrt(b);
	
    else
        
        % auxiliar variables
        i = sqrt(-1);
        b = abs(b);
        A = 0.5*sqrt(pi)*a/sqrt(b);
        xminus = sqrt(b)*(c-xmin);
        xplus  = sqrt(b)*(c-Xsupp);
        
        % CDF function
        Xcdf = A*(erfi(xminus)-erfi(xplus));
        
        % quantile function
        %Xcdfinv = c + i*erfinv(i*(-Xsupp/A + erfi(xminus)))/sqrt(b);
    end
    
    % quantile function support
    Xprob = linspace(0,1,Nx)';

    % quantile function
    Xcdfinv = interp1(Xcdf,Xsupp,Xprob,'linear','extrap')';
    
    % Entropy
    Entropy = - trapz(Xsupp,Xpdf.*log(abs(eps+Xpdf)));
    
    % PDF area
    Area = trapz(Xsupp,Xpdf);


return
% -----------------------------------------------------------------