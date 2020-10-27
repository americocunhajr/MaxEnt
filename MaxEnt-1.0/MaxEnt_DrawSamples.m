
% -----------------------------------------------------------------
%  MaxEnt_DrawSamples.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% ----------------------------------------------------------------- 
%  This functions (numerically) computes samples from a given 
%  distribution via the inverse transform method.
%
%  input:
%  Xsupp - (Nx x 1) MaxEnt PDF supoort
%  Xcdf  - (Nx x 1) MaxEnt CDF
%  Ns    - number of samples
%
%  output:
%  Xsamp - (Ns x 1) MaxEnt distribution samples
% ----------------------------------------------------------------- 

% -----------------------------------------------------------------
function Xsamp = MaxEnt_DrawSamples(Xsupp,Xcdf,Ns)

    % check number of arguments
    if nargin < 3
        error('Too few inputs.')
    elseif nargin > 3
        error('Too many inputs.')
    end
    
    % check arguments
    if length(Xsupp) ~= length(Xcdf)
        error('Xsupp and Xcdf vectors must be same length')
    end
    
    if Ns < 2
        error('Nx must be an integer greater than one.')
    end
    
    % uniform samples in [0,1]
    U = rand(Ns,1);
    
	% generate samples via inverse transform method
    Xsamp = interp1(Xcdf,Xsupp,U,'linear','extrap');

return
% -----------------------------------------------------------------