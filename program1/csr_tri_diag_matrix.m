function [A] = csr_tri_diag_matrix(n)
    % csr_tri_diag_matrix - Create a triangle diagonal Compressed Sparse Row matrix with specific requirements given in the problem.
    % Description:
    %     Modified that the first non-zero element of each line is set to be diagonal element.
    % 
    % Inputs:
    %     n - the size of matrix a
    % 
    % Outputs:
    %     A - the csr matrix structure with:
    %         A.val - column vector storing non-zero values
    %         A.col_ind - column vector storing column index in A.val
    %         A.row_ptr - column vector storing the first non-zero element in each line,
    %             and the last (n+1) element is num of non-zero elements + 1.

    %% initialize A to preallocate memory
    NUM_VAL = 3 * n - 2;
    A.val = zeros(NUM_VAL, 1);
    A.col_ind = zeros(NUM_VAL, 1);
    A.row_ptr = zeros(n + 1, 1);
    step_num_val = 0; % the number of implemented elements of A.val for each iteration

    %% fill in matrix A
    for row = 1 : n
        A.row_ptr(row) = step_num_val + 1;

        % set main diagonal element
        step_num_val = step_num_val + 1;
        A.val(step_num_val) = 2; 
        A.col_ind(step_num_val) = row;

        % set lower diagonal element
        if row ~= 1
            step_num_val = step_num_val + 1;
            A.val(step_num_val) = -1;
            A.col_ind(step_num_val) = row - 1;
        end

        % set upper diagonal element
        if row ~= n
            step_num_val = step_num_val + 1;
            A.val(step_num_val) = -1;
            A.col_ind(step_num_val) = row + 1;
        end
    end
    A.row_ptr(n + 1) = step_num_val + 1;

    A.val = A.val .* n;
end
