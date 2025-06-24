function jacob_func = dBdu(Nx, Nu, dfdu)
    import casadi.*

    % 定义符号变量
    X = casadi.MX.sym('X', Nx, 1);
    U = casadi.MX.sym('U', Nu, 1);

    % 调用 MEX 编译的函数，返回的是一个 MX 对象
    fu = dfdu(X, U);
    fu = reshape(fu,numel(fu),1);

    % 数值差分计算雅可比矩阵
    epsilon = 1e-3;
    jacob = MX.zeros(Nx*Nu, Nu);  % 预分配符号矩阵


    for i = 1:Nu
        % 对 X 进行扰动
        AA = zeros(Nu,1);
        AA(i) = epsilon;
        
        % 计算扰动后的 f_x_pre

        f_u_perturbed = dfdu(X ,U+AA);
        f_u_perturbed = reshape(f_u_perturbed,numel(f_u_perturbed),1);
        
        % 计算数值差分
        jacob(:, i) = (f_u_perturbed - fu) / epsilon;  % 使用数值差分计算雅可比矩阵
    end

    % 转换为 MEX 可调用的函数
    jacob_func = Function('dBdu', {X, U}, {jacob});
end