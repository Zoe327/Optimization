%Finite Difference Algorithm 
function [ grad_V ] = finiteDiff( f,v, xj)

% initalize the value of step size h with default value of 10^-7
N = length(v);
h = 10^-6;
grad_V = zeros(1, N);

% loop through every single element
for i = 1:N
    % create the increment of xj
    dx = zeros(N,1);
    dx(i,1) = h;
    xj_new = xj + dx;
    % calculate the finite difference
    grad_V(1,i) = (subs(f,v,xj_new) - subs(f,v,xj))./h;
end

end