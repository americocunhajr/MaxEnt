
% -----------------------------------------------------------------
%  MaxEnt_Verification2_TruncExp2.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% ----------------------------------------------------------------- 
%  This script executes verification tests to check the accucary
%  of a numerical scheme implementation to compute the MaxEnt
%  distribution for the case where the following set of independet
%  statistical information is provied:
%  - finite support: [xmin,xmax]
%  - mean value: mu1
%  
%  Remark 1:
%  The MaxEnt distribution is a 2-parameters truncated exponential
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

% check for consistency
if mu1 <= xmin || mu1 >= xmax
    error('mu1 must be in (xmin,xmax) interval');
end
% -----------------------------------------------------------------


% compute MaxEnt distribution
% -----------------------------------------------------------------

% statistical moments values vector
b = [1; mu1];

% compute MaxEnt distribution 
[lambda,Xpdf,Xsupp,Xcdf,Xcdfinv,Xprob,Entropy,Area] = ...
                          MaxEnt_MomConstr(xmin,xmax,Nx,b);

% reference distribution
[Xpdf_ref,Xsupp_ref,Xcdf_ref,Xcdfinv_ref,Xprob_ref,Entropy_ref,Area_ref] = ...
                                    MaxEnt_TruncExp2(xmin,xmax,lambda,Nx);
% -----------------------------------------------------------------


% post-processing
% -----------------------------------------------------------------

% known information
disp(' ')
disp('Known information:')
disp(['* support = [',num2str(xmin),' ',num2str(xmax),']'])
disp(['* mean    = ' ,num2str(mu1)                       ])

% MaxEnt distribution
disp(' ')
disp('MaxEnt Dist:')
disp('Truncated Exponential (2 parameters)')
                      
% Lagrange multipliers
disp(' ')
disp('Lagrange multipliers:')
disp(lambda)

% entropy difference
disp(' ')
disp('entropy difference:')
disp(abs(Entropy-Entropy_ref))

% area difference
disp(' ')
disp('area difference:')
disp(abs(Area-Area_ref))

% plot PDF curve
figure(1)
plot(Xsupp,Xpdf,'b','LineWidth',3)
hold on
plot(Xsupp_ref,Xpdf_ref,'--r','LineWidth',3)
hold off
xlabel('support')
ylabel('probability density function')
title('Verification Test')
legend('MaxEnt code','MATLAB','Location','best')
xlim([xmin xmax])

% plot PDF curve
figure(2)
plot(Xsupp,Xcdf,'b','LineWidth',3)
hold on
plot(Xsupp_ref,Xcdf_ref,'--r','LineWidth',3)
hold off
xlabel('support')
ylabel('cumulative distribution function')
title('Verification Test')
legend('MaxEnt code','MATLAB','Location','best')
xlim([xmin xmax])
ylim([0 1])

% plot quantile function curve
figure(3)
plot(Xprob,Xcdfinv,'b','LineWidth',3)
hold on
plot(Xprob_ref,Xcdfinv_ref,'--r','LineWidth',3)
hold off
xlabel('support')
ylabel('quantile function')
title('Verification Test')
legend('MaxEnt code','MATLAB','Location','best')
xlim([0 1])
ylim([xmin xmax])
% -----------------------------------------------------------------