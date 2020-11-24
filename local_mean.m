function a = local_mean(x,lambda)
% The function "local_mean.m" provides an implementation of the local mean 
% estimation introduced in: MA Colominas, G Schlotthauer, ME Torres, 
% "An unconstrained optimization approach to
% Empirical Mode Decomposition", Digital Signal Processing 40 (2015) 164-175. 
% 
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
% Author: Marcelo A. Colominas (email: macolominas@bioingenieria.edu.ar)
% March 2015
    x = x(:);
    e = ones(length(x),1);
    D = spdiags([e -2*e e], -1:1, length(x), length(x));
    D2 = D'*D;
    extremos = find(diff(x(1:end-1)).*diff(x(2:end))<0) + 1;
        P = spalloc(length(x),length(x),3*length(extremos));
    for j = 2:length(extremos)-1
        P(extremos(j),extremos(j-1)) = (extremos(j+1)-extremos(j)) / (extremos(j+1) - extremos(j-1));
        P(extremos(j),extremos(j)) = 1;
        P(extremos(j),extremos(j+1)) = (extremos(j)-extremos(j-1)) / (extremos(j+1) - extremos(j-1));
     end;
    P2 = P'*P;
    a = (P2 + lambda*D2)\(P2*x);
    