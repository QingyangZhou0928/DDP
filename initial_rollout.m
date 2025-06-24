function [xtraj, utraj, J] = initial_rollout(Nt, xtraj,utraj, xgoal, Q, QN,R,dynamics_rk4_step)
%% Initial Rollout
    import casadi.*
    %xtraj = MX(xtraj);
    for k = 1:(Nt-1)
        xtraj(:,k+1) = full(dynamics_rk4_step(xtraj(:,k), utraj(:,k)));  % 一步RK4
    end
    J = cost(xtraj, utraj, Nt, xgoal, Q, QN,R);
end