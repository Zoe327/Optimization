function [ V_min, x_opt, Iter, Cxeq, xj_set] = augLagrangian( f,v...
    ,x0, Ceq, Cineq, threshold, lambda, mu, beta, gamma )

xj = x0;
run = true;
j = 1;
ithreshold = 10^-3;
compZero = zeros(size(Ceq,1),1);
n = size(Cineq,2);
xj_set = xj;

while run
    % calculate the value of ineq at xj
    ModCineq = Cineq;
    tempC = 0.5.*beta + gamma.*subs(Cineq, v, xj);
    
    % if the g(x) <= 0, then the ineq is met, remove that constraint
    for i = 1:n
        if (tempC(i,1) <= 0)
            ModCineq(i,1) = 0;
        end
    end
    Cineqval = subs(Cineq, v, xj)
    % construct the augmented Largrangian function
    Phi = f + (mu)*sum(Ceq.^2)/2 - lambda.'*Ceq ...
        + sum(ModCineq.^2,1)
    % minimize the function using SD
    [Vi_min, Xi_opt, TempIter, v_set, Xji_set] =...
        CG(Phi, v, xj, 0, ithreshold);
    
    % update the point position
    xj = Xi_opt;
    xj_set = [xj_set xj];
 
    % compute the constraint
    CeqValue = subs(Ceq,v,xj)
    
    if (norm(CeqValue) < threshold)
        run = false;
    else
        % update the parameters
        beta = beta + 2*gamma.*subs(Cineq, v, xj);
        beta = max(beta, compZero);
        gamma = gamma*2^j;
        lambda = lambda - mu*CeqValue;
        j = j+1;
        mu = mu*(2^j);
    end
    
end

V_min = subs(f,v,xj);
x_opt = xj;
Iter = j;
Cxeq = CeqValue;

end
