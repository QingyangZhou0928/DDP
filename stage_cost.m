function stagecost = stage_cost(x, u, xgoal, Q, R)
% define stage cost
    stagecost = 0.5*((x-xgoal)'*Q*(x-xgoal)) + 0.5*u'*R*u;
end