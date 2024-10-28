function dst = csr_gd_iteration(A, b, src)
    % Iterate A * dst(k + 1) = b, src = dst(k) using gradient descent iteration.
    % Called function "csr_vmult()".
    % 
    % Description:
    %     gradient descent iteration for Ax=b:
    %     $r$ is negative gradient direction
    %     $\alpha$ is step size
    %     $(,)$ is vector norm
    %     r(k) = b - A * x(k)
    %     \alpha(k) = \frac{(r(k), r(k))}{(A * r(k), r(k))}
    %     x(k + 1) = x(k) + \alpha(k) * r(k)
    % 
    %     So dst = x

    r = b - csr_vmult(A, src);
    alpha = (r' * r) / (r' * csr_vmult(A, r));
    dst = src + alpha * r;

end
