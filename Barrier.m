function [ V_min, x_opt, Iter, xj_set ]...
    = Barrier( f,v, Cineq, k, x0, threshold )
j = 1;
xj = x0;
run = true;
threshold = 10^-5;
n = size(Cineq,1);
xj_set = xj;

while run
    r = 1/(k)^j;
    % convert constraint problem to an unconstraint one
    Phi = f - r*sum(Cineq.^-1,1);
    % solve for the new unconstrainted problem using secant method
    [new_xj, TempIter] = Secant(Phi, v, xj, 1, 10^-4);
    % As barrier method is an interior point method and the cost
    % function blows up if it gets close to the boundary.
    if (norm(new_xj - xj) < threshold)
        run = false;
    else
        xj = new_xj;
        xj_set = [xj_set xj];
        j = j+1;
    end
end
V_min = subs(f,v,new_xj);
x_opt = new_xj;
Iter = j;
end