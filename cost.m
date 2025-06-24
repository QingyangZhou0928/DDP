function J = cost(xtraj, utraj, Nt, xgoal, Q, QN,R)
    J = 0;
    for k = 1:(Nt-1)
        J = J+stage_cost(xtraj(:,k),utraj(:,k), xgoal, Q, R);
    end
    J = J+terminal_cost(xtraj(:,Nt), xgoal, QN);
end
