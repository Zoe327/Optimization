
%% Unconstrained Problems Testing 
syms x1 x2 x3 x4 x5 x6;
v = [x1; x2; x3; x4; x5; x6];
C = [9 1 7 5 4 7; 1 11 4 2 7 5; 7 4 13 5 0 7; 5 2 5 17 1 9; ...
    4 7 0 1 21 15; 7 5 7 9 15 27];
f = 5+ [1 4 5 4 2 1]*v + v.'*C*v;

% initalize the starting point
x0 = [1 1 0 0 1 1]';
threshold = 10^-4;
% Steepest Decent 
% jacobian
[Voptm1, Xoptm1, Iter1, v_set, xj_set] = a_SD(f,v,x0,0, threshold);
% finite difference 
[Voptm2, Xoptm2, Iter2, v_set, xj_set] = a_SD(f,v, x0, 1, threshold);

% Conjugate Gradient 
 [Voptm3, Xoptm3, iter3, v_set, xj_set] = CG(f, v, x0, 1, threshold);
  
%Secant Algorithm
% jacobian
 [Voptm5, Xoptm5, Iter5, v_set, xj_set] = Secant(f,v,x0,0, threshold);
% finite difference method
[Voptm6, Xoptm6, Iter6, v_set, xj_set] = Secant_debug(f,v,x0,1,threshold);


%% Constraint problems testing
syms x1 x2
x = [x1; x2];
V = abs(x1-2) + abs(x2-2);

Eq = x1^2+ x2^2-1;
Ineq =  x1- x2^2;
x0 = [1.5, 1.5]';
sigma = 10;
stop_param = 10^-4;

% penalty method
tic
 [V_min, X_min, iteration, xj_set] = ...
    Penalty_mix(V,x,sigma,Eq,Ineq,x0,stop_param);
toc
V_min_penalty = double(V_min);
X_min_penalty = double(X_min);


%Augment lagrangian
lambda = 10;
mu = 10;
beta = 0;
gamma = 2;

tic
[V_min, X_min, iteration, xj_set] = AugLagrangian(V,x,x0,Eq,...
     Ineq,stop_param,lambda,mu,beta,gamma);
toc
V_min_auglag = double(V_min);
X_min_auglag = double(X_min);
