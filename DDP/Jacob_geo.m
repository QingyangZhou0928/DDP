function J = Jacob_geo(q,DH)
% 输入：各转动关节角度参数q(3*1)
% 输出：末端雅可比阵J
    import casadi.*
    n = length(q);
    T = dh(q,DH);
    J = MX(zeros(6, n));
    J(:, 1) = [cross([0; 0; 1], T(1:3, 4*n)); [0; 0; 1]];
    for i = 2:n
        J(:, i) = [cross(T(1:3,4*i-5), (T(1:3, 4*n)-T(1:3, 4*i-4))); T(1:3, 4*i-5)];
    end
end