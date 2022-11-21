%% Compute P_t given Q for heterogeneous insertion rates

% Here, we compute the transition probability matrix P given the rate 
% matrix Q. This rate matrix models the accumulation of 5 edits (of any type).
% It assumes that the rates r1, ..., r5 are all pairwise different. If not,
% the transition rate matrix Q is not of full rank and cannot be 
% diagonalised, i.e. the transition probability matrix cannot be computed.

syms r1 r2 r3 r4 r5 t;
format long;
%r1=0.2;r2=0.1;r3=0.11;r4=0.105;r5=0.8;t=6;

Q = [   -r1 r1 0 0 0 0;
               0 -r2 r2 0 0 0;
               0 0 -r3 r3 0 0;
               0 0 0 -r4 r4 0;
               0 0 0 0 -r5 r5;
               0 0 0   0  0 0];
           
[V, D] = eig(Q);   
D = diag(exp(eig(Q) * t));
P = V * D * V^-1
