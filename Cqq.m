function C = Cqq(q, qdot,DH,m,I)
    import casadi.*
    n = length(q);
    h = 1e-6;
    C = MX(zeros(n, n));
    
    for i = 1:n
        for j = 1:n
            for k = 1:n
                % 计算∂M_ij/∂q_k
                theta_plus = q;
                theta_plus(k) = theta_plus(k) + h;
                M_plus = Mq(theta_plus,m,I,DH);
                
                theta_minus = q;
                theta_minus(k) = theta_minus(k) - h;
                M_minus = Mq(theta_minus,m,I,DH);
                
                dM_ij = (M_plus(i,j) - M_minus(i,j))/(2*h);
                
                % 计算∂M_jk/∂q_i
                theta_plus = q;
                theta_plus(i) = theta_plus(i) + h;
                M_plus = Mq(theta_plus,m,I,DH);
                
                theta_minus = q;
                theta_minus(i) = theta_minus(i) - h;
                M_minus = Mq(theta_minus,m,I,DH);
                
                dM_jk = (M_plus(j,k) - M_minus(j,k))/(2*h);
                
                C(i,j) = C(i,j) + (dM_ij - 1/2 * dM_jk) * qdot(k);
            end
        end
    end
end