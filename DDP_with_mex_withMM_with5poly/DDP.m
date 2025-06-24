function [xtraj,utraj,J,Jhist] = DDP(Nx, Nu, Nt, xtraj,utraj, xgoal, Q, QN,R)
% basic DDP
    import casadi.*

%% initial rollout
    J = cost(xtraj, utraj, Nt, xgoal, Q, QN,R);
%% initial parameters
    d = ones(Nu, Nt-1);
    iter = 0;

%% 大循环
    while max(abs(d(:))) > 1e-3
        iter = iter+1;
        disp(iter)
    % Backward Pass
        [delta_J,d,K] = backward_pass_DDP(Nx,Nt,Nu, xtraj, xgoal,utraj,QN,Q,R,iter);
    % Forward rollout with line search
        xn(:,1) = xtraj(:,1);
        alpha = 1.0;
        while true
            for k = 1:(Nt-1)
                un(:,k) = utraj(:,k) - alpha*d(:,k) - K(:,:,k)*(xn(:,k)-xtraj(:,k));
                xn(:,k+1)= full(dynamics_rk4_step(xn(:,k),un(:,k)));
            end
            Jn = cost(xn, un, Nt, xgoal, Q, QN,R);
    
            if isnan(Jn) || Jn > (J - 1e-2*alpha*delta_J)
                alpha = 0.5*alpha;
            else
                break
            end
        end
        J = Jn;
        xtraj = xn;
        utraj = un;
        Jhist(iter) = J;
    end
    disp("DDP法迭代次数：");
    disp(iter);
end