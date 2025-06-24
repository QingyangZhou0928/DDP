function dynamics_rk4_step = rk4_step( h,dynamics_jiexi,Nu)
    import casadi.*
    % 一步RK4积分
    % 定义符号变量
    x = MX.sym('x', 2*Nu+3, 1);  % 状态变量
    u = MX.sym('u', Nu, 1);  % 控制输入

    k1 = dynamics_jiexi(x, u);
    k2 = dynamics_jiexi(x + 0.5*h*k1, u);
    k3 = dynamics_jiexi(x + 0.5*h*k2, u);
    k4 = dynamics_jiexi(x + h*k3, u);

    x_next = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    dynamics_rk4_step = Function('dynamics_rk4_step', {x, u}, {x_next});
end
