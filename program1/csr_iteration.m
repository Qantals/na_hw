function dst = csr_iteration(A, b, src, method)
    % csr_iteration - Iteration method cluster for one step.
    % 
    % Description:
    % implemented for some methods, seen subfunctions below in this function file.
    % 
    % Inputs:
    %     A - cofficient matrix in Ax=b.
    %     b - result column vector in Ax=b.
    %     src - last step like x(k).
    %     method - choices between "Jacobi", "Gauss-Seidel" and "gradient descent"
    % 
    % Outputs:
    %     dst - solution like x(k + 1)

    method = string(method);
    switch method
        case "Jacobi"
            dst = csr_jacobi_iteration(A, b, src);
        case "Gauss-Seidel"
            dst = csr_gs_iteration(A, b, src);
        case "gradient descent"
            dst = csr_gd_iteration(A, b, src);
        otherwise
            error("Invalid method. Choose 'Jacobi', 'Gauss-Seidel' or 'gradient descent'.");
    end
end


function dst = csr_jacobi_iteration(A, b, src)
    % csr_jacobi_iteration - Iterate A * dst(k + 1) = b, src = dst(k) using Jacobi iteration.
    % 
    % Description:
    %     Jacobi iteration for Ax=b:
    %     A = D - L - U
    %     x(k + 1) = D^{-1} * (L + U) * x(k) + D^{-1} * b
    % 
    %     From code to result above:
    %     dst(k + 1) = (b + (L + U) * src) * D^{-1}
    %          = b * D^{-1} + (L + U) * src * D^{-1}
    %          = D^{-1} * b + D^{-1} * (L + U) * src
    %          = D^{-1} * (L + U) * src + D^{-1} * b 
    % 
    %     So dst(k) = src

    n = length(src);
    dst = b;

    for i = 1 : n
        j = A.row_ptr(i) + 1 : A.row_ptr(i + 1) - 1; % note: this is L and U without D!
        dst(i) = (dst(i) - A.val(j)' * src(A.col_ind(j))) / A.val(A.row_ptr(i));
    end

end

function dst = csr_gs_iteration(A, b, src)
    % csr_gs_iteration - Iterate A * dst(k + 1) = b, src = dst(k) using Gauss-Seidel iteration.
    % 
    % Description:
    %     Gauss-Seidel iteration for Ax=b:
    %     A = D - L - U
    %     x(k + 1) = (D - L)^{-1} * U * x(k) + (D - L)^{-1} * b
    % 
    %     From result above to code:
    %     dst(k + 1) = (D - L)^{-1} * U * dst(k) + (D - L)^{-1} * b
    %                = (D - L)^{-1} * (b + U * dst(k))
    %     We can left multiply (D - L) and use lower triangle characteristic to do iteration.

    n = length(src);
    dst = b;

    for i = 1 : n
        j = A.row_ptr(i) : A.row_ptr(i + 1) - 1; % including D, L, U
        j_L = j(A.col_ind(j) < i);
        j_U = j(A.col_ind(j) > i);
        j_D = A.row_ptr(i);
        dst(i) = ((dst(i) - A.val(j_U)' * src(A.col_ind(j_U)))...
                    - dst(A.col_ind(j_L))' * A.val(j_L))...
                    / A.val(j_D);
    end

end

function dst = csr_gd_iteration(A, b, src)
    % csr_gd_iteration - Iterate A * dst(k + 1) = b, src = dst(k) using gradient descent iteration.
    %     Called function "csr_vmult()".
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
