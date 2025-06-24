function M = Mq(q,m,I,DH)
% 输入：各转动关节角度参数q(n*1)(若n表示关节数量的话)，连杆质量m(n*1)，转动惯量矩阵I(n*1)
% 输出：拉格朗日方程中二阶导项系数矩阵B(q)
    import casadi.*
    n = length(q);
    M = zeros(n, n);
    J = Jacob_geo(q,DH);
    for i = 1:n
        Ji = MX(zeros(6, n));
        for j = 1:i
            Ji(:, j) = J(:, j);
        end
        M = M+m(i)*Ji(1:3, :)'*Ji(1:3, :)+Ji(4:6, :)'*I(3*i-2:3*i,:)*Ji(4:6, :);
    end
end