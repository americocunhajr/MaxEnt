
% -----------------------------------------------------------------
%  MaxEnt_Example5_Beta.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% ----------------------------------------------------------------- 
%  This script ilustrates how to numerically compute the MaxEnt
%  distribution for the case of where the following set of 
%  independet statistical information is provied:
%  - finite support: [0,1]
%  - geometric mean at the support lower bound: mu_log_0
%  - geometric mean at the support upper bound: mu_log_1
%  
%  Remark 1:
%  The MaxEnt distribution is the standard beta
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
xmax = 1.0;

% mean (xmin <= mu1 <= xmax)
mu1 = 0.55;

% coeficient of variation (cv > 0)
cv = 0.25;

% standard deviation (sigma > 0)
sigma = cv*mu1;

% shape parameters adition
nu = mu1*(1-mu1)/(sigma^2)-1;

% shape parameter 1
shape1 = mu1*nu;

% shape parameter 2
shape2 = (1-mu1)*nu;

% geometric mean
mu_log_0 = psi(shape1)-psi(shape1+shape2);
mu_log_1 = psi(shape2)-psi(shape1+shape2);

% check for consistency
if mu1 <= xmin || mu1 >= xmax
    error('mu1 must be in (xmin,xmax) interval');
end
% -----------------------------------------------------------------


% compute MaxEnt distribution
% -----------------------------------------------------------------

% statistical moments values vector
b = [1; mu_log_0; mu_log_1];

% statistical properties values vector
gfunc = @(x) StatPropFunc(x,xmin,xmax);

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
disp(['* support        = [',num2str(xmin),' ',num2str(xmax),']'])
disp(['* mean           = ' ,num2str(mu1)                       ])
disp(['* geo. mean at 0 = ' ,num2str(mu_log_0)                  ])
disp(['* geo. mean at 1 = ' ,num2str(mu_log_1)                  ])

% MaxEnt distribution
disp(' ')
disp('MaxEnt Dist:')
disp('Beta')

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
title('MaxEnt distribution: Beta')
xlim([xmin xmax])

% plot CDF curve
figure(2)
plot(Xsupp,Xcdf,'r','LineWidth',3)
xlabel('support')
ylabel('cumulative distribution function')
title('MaxEnt distribution: Beta')
xlim([xmin xmax])
ylim([0 1])

% plot CDFinv curve
figure(3)
plot(Xprob,Xcdfinv,'g','LineWidth',3)
xlabel('probability')
ylabel('quantile function')
title('MaxEnt distribution: Beta')
xlim([0 1])
ylim([xmin xmax])

% plot samples
figure(4)
plot(1:Ns,Xsamp,'xb');
xlabel('sample index')
ylabel('sample value')
title('MaxEnt distribution: Beta')
xlim([1 Ns])
% -----------------------------------------------------------------

% statistical properties function
% -----------------------------------------------------------------
function g = StatPropFunc(x,xmin,xmax)
    % mesh size
    Nx = length(x);
    % preallocate memory for the constraint functions
    g = zeros(Nx,3);
    % constraint functions
    g(:,1) = ones(Nx,1);
    g(:,2) = log(abs((x-xmin)+eps));
    g(:,3) = log(abs((xmax-x)-eps));
end
% -----------------------------------------------------------------