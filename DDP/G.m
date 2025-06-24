function g = G(q,g0,m,DH)
% 输入：各转动关节角度参数q(n*1)，连杆质量m(n*1)
% 输出：拉格朗日方程中常数项g(q)
    import casadi.*
    n = length(q);
    g = MX(zeros(n, 1));
    J = Jacob_geo(q,DH);
    for i = 1:n
        for j = 1:n
            if ((j-1)>=i)
                g(i) = g(i)-0.5*m(i)*g0'*J(1:3, i);
            end
            if (j>=i)
                g(i) = g(i)-0.5*m(i)*g0'*J(1:3, i);
            end
        end
    end
end
