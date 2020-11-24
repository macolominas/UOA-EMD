# UOA-EMD
A MATLAB file implementing the adaptive local mean proposed in: MA Colominas, G Schlotthauer, ME Torres. "An Unconstrained Optimization Approach to Empirical Mode Decomposition", Digital Signal Processing (Elsevier), vol. 40, pp. 164-175, 2015.

% The function "local_mean.m" provides an implementation of the local mean 
% estimation introduced in: MA Colominas, G Schlotthauer, ME Torres, 
% "An unconstrained optimization approach to
% Empirical Mode Decomposition", Digital Signal Processing 40 (2015) 164-175. 
% 
% It implements Eqs. (5) from Algorithm 1. Following this algorithm, the whole method
% can be easily implemented.
%
% The optimization problem to be solved is:
%
%  min_a ||P (x-a)||^2 + lambda * ||D a||^2,
%
% where x is the original signal, a its estimated local mean, P is a matrix
% constructed from the local extrema of x, D is a second-order difference
% matrix, ||.|| stands for the Euclidian norm, and lambda>0 is a
% regularization parameter (recommended value: lambda = 1). The problem is
% designed so the signal x-a approximately fulfills IMF conditions.
%
% SYNTAX a = local_mean(x,lambda)
% 
% INPUT: x (signal)
%        lambda (regularization parameter)
%
% OUTPUT: a (estimated local mean)
%
% Version 1.0 for Matlab R2013b
% Author: Marcelo A. Colominas (email: macolominas@conicet.gov.ar)
% March 2015
