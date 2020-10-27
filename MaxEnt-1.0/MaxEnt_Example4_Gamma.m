
% -----------------------------------------------------------------
%  MaxEnt_Example4_Gamma.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% ----------------------------------------------------------------- 
%  This script ilustrates how to numerically compute the MaxEnt
%  distribution for the case of where the following set of 
%  independet statistical information is provied:
%  - infinite positive support: [0,infty]
%  - mean value: mu1
%  - geometric mean at the support lower bound: mu_log_0
%  
%  Remark 1:
%  The MaxEnt distribution is the gamma
%  
%  Remark 2:
%  This code can be easly adapted to obtain the MaxEnt distribution
%  associated with other types of statistical information provided
%  by the user.
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
% -----------------------------------------------------------------

clc
clear
close all


% program header
% -----------------------------------------------------------
disp('================================================')
disp('   MaxEnt - Maximum Entropy Code                ')
disp('   by A. Cunha Jr                               ')
disp('                                                ')
disp('   This is an easy to run code for calculation  ')
disp('   of the MaxEnt distribution given a set of    ')
disp('   independent set of statistical information.  ')
disp('================================================')
% -----------------------------------------------------------


% stochastic simulation paramters
% -----------------------------------------------------------------

% number of points for support discretization (Nx > 1)
Nx = 1000;

% number of random samples
Ns = 2048;

% number of bins
Nbins = round(sqrt(Ns));
% -----------------------------------------------------------------


% known information
% -----------------------------------------------------------------

% support lower bound
xmin = 0.0;

% support upper bound
xmax = 2.0;

% mean (xmin <= mu1 <= xmax)
mu1 = 0.3;

% coeficient of variation (0 < cv < 1/sqrt(2))
cv = 0.65;

% standard deviation (sigma > 0)
sigma = cv*mu1;

% shape parameter
shape = (mu1/sigma)^2;

% scale parameter
scale = sigma^2/mu1;

% geometric mean
mu_log_0 = psi(shape) + log(scale);

% check for consistency
if mu1 <= xmin || mu1 >= xmax
    error('mu1 must be in (xmin,xmax) interval');
end

% check for consistency
if cv <= 0.0 
    error('cv cannot be negative or zero');
end
% -----------------------------------------------------------------


% compute MaxEnt distribution
% -----------------------------------------------------------------

% statistical moments values vector
b = [1; mu1; mu_log_0];

% statistical properties values vector
gfunc = @(x) StatPropFunc(x);

% compute MaxEnt distribution 
[lambda,Xpdf,Xsupp,Xcdf,Xcdfinv,Xprob,Entropy,Area] = ...
                          MaxEnt_GenConstr(xmin,xmax,Nx,b,gfunc);
    
% generate samples (inverse transform method)
Xsamp = MaxEnt_DrawSamples(Xsupp,Xcdf,Ns);

% histogram
[Xbins,Xfreq] = randvar_pdf(Xsamp,Nbins);
% -----------------------------------------------------------------


% post-processing
% -----------------------------------------------------------------

% known information
disp(' ')
disp('Known information:')
disp(['* support   = [',num2str(xmin),' ',num2str(xmax),']'])
disp(['* mean      = ' ,num2str(mu1)                       ])
disp(['* geo. mean = ' ,num2str(mu_log_0)                  ])

% MaxEnt distribution
disp(' ')
disp('MaxEnt Dist:')
disp('Gamma')
                      
% Lagrange multipliers
disp(' ')
disp('Lagrange multipliers:')
disp(lambda)

% plot PDF curve
figure(1)
bar(Xbins,Xfreq,'m');
hold on
plot(Xsupp,Xpdf,'b','LineWidth',3)
hold off
xlabel('support')
ylabel('probability density function')
title('MaxEnt distribution: Gamma')
xlim([xmin xmax])

% plot CDF curve
figure(2)
plot(Xsupp,Xcdf,'r','LineWidth',3)
xlabel('support')
ylabel('cumulative distribution function')
title('MaxEnt distribution: Gamma')
xlim([xmin xmax])
ylim([0 1])

% plot CDFinv curve
figure(3)
plot(Xprob,Xcdfinv,'g','LineWidth',3)
xlabel('probability')
ylabel('quantile function')
title('MaxEnt distribution: Gamma')
xlim([0 1])
ylim([xmin xmax])

% plot samples
figure(4)
plot(1:Ns,Xsamp,'xb');
xlabel('sample index')
ylabel('sample value')
title('MaxEnt distribution: Gamma')
xlim([1 Ns])
% -----------------------------------------------------------------

% statistical properties function
% -----------------------------------------------------------------
function g = StatPropFunc(x)
    % mesh size
    Nx = length(x);
    % preallocate memory for the constraint functions
    g = zeros(Nx,3);
    % constraint functions
    g(:,1) = ones(Nx,1);
    g(:,2) = x;
    g(:,3) = log(abs(x+eps));
end
% -----------------------------------------------------------------