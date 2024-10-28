function [dst, r_n, p_n] = csr_cg_iteration(A, src, r, p)
    % csr_cg_iteration - Iterate A * x(k + 1) = b, src = x(k), dst = x(k + 1),
    %     r = r(k), p = p(k) using conjugate gradient iteration.
    %     Called function "csr_vmult()".
    % 
    % Description:
    %     gradient descent iteration for Ax=b:
    %     $r$ is negative gradient direction
    %     $p$ is conjugate direction
    %     $\alpha$ is step size
    %     %\beta% is cofficient of $p(k)$ to construct conjugate direction $p(k + 1)$
    %     $(,)$ is vector norm
    % 
    %     \alpha(k) = \frac{(r(k), r(k))}{(A * p(k), p(k))}
    %     x(k + 1) = x(k) + \alpha(k) * p(k)
    %     r(k + 1) = r(k) - \alpha(k) * A * p(k)
    %     \beta(k) = \frac{r(k + 1), r(k + 1)}{r(k), r(k)}
    %     p(k + 1) = r(k + 1) + \beta(k) * p(k)

    alpha = (r' * r) / (p' * csr_vmult(A, p));
    dst = src + alpha * p;
    r_n = r - alpha * csr_vmult(A, p);
    beta = (r_n' * r_n) / (r' * r);
    p_n = r_n + beta * p;

end