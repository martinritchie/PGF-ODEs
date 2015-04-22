% S:        fraction of susceptible individuals,
% I:        fraction of infected individuals,
% R:        fraction of recovered individuals,
% T:        time,
function [ S, I, R, T ] = episolve()
global alpha Tend
% ODE integration
options = odeset('AbsTol',1e-8,'RelTol',1e-8);
[T,Y] = ode45(@func,[0 Tend],alpha,options);
% I will always be the 2nd to last varibale
I = Y(:,end-1);
% R will always be the last.
R = Y(:,end);
S = 1 - I - R;
end

