function jacob_func = dAdu(Nx,Nu,dfdx)
% f: 离散系统xk+1=f_discrete(xk,u)
% x: 当前状态
% u: 控制输入
% epsilon: 微小扰动值，用于数值差分
    
% 调用 rk4_step 函数计算状态方程
    import casadi.*
    % 定义符号变量
    X = MX.sym('X', Nx, 1);  % 状态变量
    U = MX.sym('U', Nu, 1);  % 控制输入

    fx = dfdx(X, U);
    fx = reshape(fx, numel(fx), 1);

    % 计算 f 对 X 的雅可比矩阵
    jacob = jacobian(fx, U);  % 使用 CasADi 的 jacobian 函数
    
    % 将雅可比矩阵转换为数值函数，方便后续调用
    jacob_func = Function('dAdu', {X, U}, {jacob});

end

