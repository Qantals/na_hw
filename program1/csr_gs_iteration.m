function dst = csr_gs_iteration(A, b, src)
    % Iterate A * dst(k + 1) = b, src = dst(k) using Gauss-Seidel iteration.
    % 
    % Description:
    %     % Gauss-Seidel iteration for Ax=b:
    %     A = D - L - U
    %     x(k + 1) = (D - L)^{-1} * U * x(k) + (D - L)^{-1} * b
    % 
    %     % From result above to code:
    %     dst(k + 1) = (D - L)^{-1} * U * dst(k) + (D - L)^{-1} * b
    %                = (D - L)^{-1} * (b + U * dst(k))
    %     % We can left multiply (D - L) and use lower triangle characteristic to do iteration.

    n = length(src);
    dst = b;

    for i = 1 : n
        j = A.row_ptr(i) : A.row_ptr(i + 1) - 1; % including D, L, U
        j_L = j(A.col_ind(j) < i);
        j_U = j(A.col_ind(j) > i);
        j_D = A.row_ptr(i);
        % dst(i) = ((dst(i) - A.val(j_U)' * src(A.col_ind(j_U))) - dst(j_L)' * A.val(j_L)) / A.val(j_D);
        dst(i) = ((dst(i) - A.val(j_U)' * src(A.col_ind(j_U)))...
                    - dst(A.col_ind(j_L))' * A.val(j_L))...
                    / A.val(j_D);
    end

end
