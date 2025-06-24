function isPD = isPositiveDe(A)
    isPD = false; % 默认假设不是正定的
    try
        chol(A); % 尝试 Cholesky 分解
        isPD = true; % 如果成功则是正定的
    catch
        isPD = false; % 如果失败则不是正定的
    end
end