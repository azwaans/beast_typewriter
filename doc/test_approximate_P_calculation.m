
% Is approximate calculation of P accurate?

% Define rate matrix Q in terms of a single rate of edit introduction, r1.
% Let t be time.
syms r1 t;

Q = [   -r1 r1 0 0;
               0 -r1 r1 0;
               0 0 -r1 r1;
               0 0 0   0];
           
% Then, we calculate the transition probability matrix, P, by approximately
% by a method in matlab
P = expm(Q*t);

% Is P correct? To test, verify that d P(t)/ dt = Q * P(t)
P_diff = diff(P);
P_prod = Q * P;

% The equality does not hold, hence P(t) is not correct
isequal(P_diff, P_prod)
