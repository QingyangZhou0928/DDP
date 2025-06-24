% 齐次变换矩阵计算函数
function T = dh(q,DH)
% 输入：各转动关节角度参数theta(n*1)
% 输出：矩阵T，T(:, :, i)表示第i个坐标系的齐次变换阵
    import casadi.*
    n = length(q);
    T = zeros(4, 4*n);
    T = MX(T);
    % 构造一个新的 DH_eff，其中 theta 列加入符号变量 q
    DH_eff = [DH(:,1) + q, DH(:,2:4)];
    
    T(:, 1:4) = rotz(DH_eff(1, 1))*tranz(DH_eff(1, 2))*tranx(DH_eff(1, 3))*rotx(DH_eff(1, 4));
    for i = 2:n
        T(:, (4*i-3):(4*i)) = T(:, (4*i-7):(4*i-4)) * rotz(DH_eff(i, 1))*tranz(DH_eff(i, 2))*tranx(DH_eff(i, 3))*rotx(DH_eff(i, 4));
    end
end