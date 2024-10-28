function dst = csr_vmult(A, src)
    % Get dst = A * src, where "A" is sparse square matrix and "src" is column vector.

    n = length(src);
    dst = zeros(n, 1);

    %% step 1: trival try with 2 for loop
    % for i = 1 : n
    %     for j = A.row_ptr(i) : A.row_ptr(i + 1) - 1
    %         dst(i) = dst(i) + A.val(j) * src(A.col_ind(j));
    %     end
    % end

    %% step 2: try with 1 for loop
    for i = 1 : n
        j = A.row_ptr(i) : A.row_ptr(i + 1) - 1;
        dst(i) = A.val(j)' * src(A.col_ind(j));
    end

    %% step 3 is not known, because dealing with sparsity is the first.
end
