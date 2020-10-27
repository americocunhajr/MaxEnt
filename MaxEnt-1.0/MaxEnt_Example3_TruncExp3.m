
% -----------------------------------------------------------------
%  MaxEnt_Example3_TruncExp3.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% ----------------------------------------------------------------- 
%  This script ilustrates how to numerically compute the MaxEnt
%  distribution for the case of where the following set of 
%  independet statistical information is provied:
%  - finite support: [xmin,xmax]
%  - mean value: mu1
%  - standard deviation: sigma
%  
%  Remark 1:
%  The MaxEnt distribution is a 3-parameters truncated exponential
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
xmin = -1.0;

% support upper bound
xmax = 1.0;

% mean (xmin <= mu1 <= xmax)
mu1 = -0.3;

% standard deviation (sigma > 0)
sigma = 0.4;

% second moment
mu2_x = mu1^2 + sigma^2;

% check for consistency
if mu1 <= xmin || mu1 >= xmax
    error('mu1 must be in (xmin,xmax) interval');
end

% check for consistency
if sigma <= 0.0 
    error('sigma cannot be negative or zero');
end
% -----------------------------------------------------------------


% compute MaxEnt distribution
% -----------------------------------------------------------------

% statistical moments values vector
b = [1; mu1; mu2_x];

% compute MaxEnt distribution 
[lambda,Xpdf,Xsupp,Xcdf,Xcdfinv,Xprob,Entropy,Area] = ...
                          MaxEnt_MomConstr(xmin,xmax,Nx,b);

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
disp(['* std. dev. = ' ,num2str(sigma)                     ])

% MaxEnt distribution
disp(' ')
disp('MaxEnt Dist:')
disp('Truncated Exponential (3 parameters)')
                      
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
title('MaxEnt distribution: Truncated Exponential (3 parameters)')
xlim([xmin xmax])

% plot CDF curve
figure(2)
plot(Xsupp,Xcdf,'r','LineWidth',3)
xlabel('support')
ylabel('cumulative distribution function')
title('MaxEnt distribution: Truncated Exponential (3 parameters)')
xlim([xmin xmax])
ylim([0 1])

% plot CDFinv curve
figure(3)
plot(Xprob,Xcdfinv,'g','LineWidth',3)
xlabel('probability')
ylabel('quantile function')
title('MaxEnt distribution: Truncated Exponential (3 parameters)')
xlim([0 1])
ylim([xmin xmax])

% plot samples
figure(4)
plot(1:Ns,Xsamp,'xb');
xlabel('sample index')
ylabel('sample value')
title('MaxEnt distribution: Truncated Exponential (3 parameters)')
xlim([1 Ns])
% -----------------------------------------------------------------
