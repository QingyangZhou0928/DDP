function dynamics_jiexi = dynamics(robot,Nu,DH)
%% 针对串联机器人
    import casadi.*
    
    m = [robot.links.m].';  % 连杆质量矩阵
    I = [robot.links.I].';  % 连杆惯性质量矩阵

    % 定义符号变量
    x = MX.sym('x', 2*Nu+4, 1);  % 状态变量
    u = MX.sym('u', Nu, 1);  % 控制输入
    q = x(1:Nu);
    qdot = x(Nu+1:2*Nu);
    end_effector = x(2*Nu+1:2*Nu+3);  % end effector末端位姿变量
    omega = x(2*Nu+4);  % MM

    % 获取动力学参数
    M_sym = Mq(q,m,I,DH);
    C_sym = Cqq(q, qdot,DH,m,I);
    g_sym = G(q,robot.gravity,m,DH);
    
    JA = Jacob_geo(q,DH);  % 分析雅可比矩阵
    JA = JA([1,2,6],:);
    %while cond(JA) > 1e12 || isnan(cond(JA))  % 条件数越大表示数值不稳定
        %warning('矩阵 JA 可能是奇异或病态的，跳过这一步或做特殊处理');
        %q=q+[0;0.01;0.02];
        %JA = Jacob_geo(q,DH);  % 分析雅可比矩阵
        %JA = JA([1,2,6],:);
        %double(JA)
    %end

    % 计算各分量
    M_1 = pinv(M_sym, 'symbolicqr');
    qddot = M_1*(u - C_sym*qdot - g_sym);  % 加速度
    end_effectordot = JA * qdot;  % 末端执行器位置变化率

    t_sample = 0.01;
    JAdelta = Jacob_geo(q+qdot*t_sample,DH);
    omegadot = (sqrt(det(JAdelta*JAdelta'))-sqrt(det(JA*JA')))/t_sample;
    % 构造 xdot
    xdot = [qdot;
            qddot;
            end_effectordot;
            omegadot];

    % 转换为数值函数，方便后续调用
    dynamics_jiexi = Function('dynamics_jiexi', {x, u}, {xdot});
end
