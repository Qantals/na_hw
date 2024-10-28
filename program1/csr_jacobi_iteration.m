function dst = csr_jacobi_iteration(A, b, src)
    % Iterate A * dst(k + 1) = b, src = dst(k) using Jacobi iteration.
    % 
    % Description:
    %     % Jacobi iteration for Ax=b:
    %     A = D - L - U
    %     x(k + 1) = D^{-1} * (L + U) * x(k) + D^{-1} * b
    % 
    %     % From code to result above:
    %     dst(k + 1) = (b + (L + U) * src) * D^{-1}
    %          = b * D^{-1} + (L + U) * src * D^{-1}
    %          = D^{-1} * b + D^{-1} * (L + U) * src
    %          = D^{-1} * (L + U) * src + D^{-1} * b 
    % 
    %     % So dst(k) = src

    n = length(src);
    dst = b;

    for i = 1 : n
        j = A.row_ptr(i) + 1 : A.row_ptr(i + 1) - 1; % note: this is L and U without D!
        dst(i) = (dst(i) - A.val(j)' * src(A.col_ind(j))) / A.val(A.row_ptr(i));
    end

end
