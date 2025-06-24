function K = comm(n, m)
% COMMUTATION MATRIX K such that vec(A') = K * vec(A) for A ∈ ℝ^{n×m}
% Input: n (rows), m (cols)
% Output: K ∈ ℝ^{nm × nm}

    K = zeros(n*m, n*m);
    for i = 1:n
        for j = 1:m
            % vec(A): A(j,i) → index (i-1)*m + j
            % vec(A'): A'(i,j) = A(j,i) → index (j-1)*n + i
            row = (j-1)*n + i;
            col = (i-1)*m + j;
            K(row, col) = 1;
        end
    end
end
