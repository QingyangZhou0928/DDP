function jacob_func = dAdx(Nx, Nu, dfdx)
    import casadi.*

    % 定义符号变量
    X = casadi.MX.sym('X', Nx, 1);
    U = casadi.MX.sym('U', Nu, 1);

    % 调用 MEX 编译的函数，返回的是一个 MX 对象
    fx = dfdx(X, U);
    fx = reshape(fx,numel(fx),1);

    % 数值差分计算雅可比矩阵
    epsilon = 1e-3;
    jacob = MX.zeros(Nx*Nx, Nx);  % 预分配符号矩阵


    for i = 1:Nx
        % 对 X 进行扰动
        AA = zeros(Nx,1);
        AA(i) = epsilon;
        
        % 计算扰动后的 f_x_pre

        f_x_perturbed = dfdx(X+AA ,U);
        f_x_perturbed = reshape(f_x_perturbed,numel(f_x_perturbed),1);
        
        % 计算数值差分
        jacob(:, i) = (f_x_perturbed - fx) / epsilon;  % 使用数值差分计算雅可比矩阵
    end

    % 转换为 MEX 可调用的函数
    jacob_func = Function('dAdx', {X, U}, {jacob});
end