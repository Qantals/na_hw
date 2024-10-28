clear; clc;
%% TASK 2
fprintf('Task 2: csr_vmult(A, b)\n');
n = 16;
h = 1 / n;
A = csr_tri_diag_matrix(n);
b = h * ones(n, 1);
csr_vmult(A, b)

%% TASK 3
fprintf('Task 3: csr_jacobi_iteration(A, b, b)\n');
csr_iteration(A, b, b, "Jacobi")

%% TASK 4
fprintf('Task 4: csr_gs_iteration(A, b, b)\n');
csr_iteration(A, b, b, "Gauss-Seidel")

%% TASK 5
fprintf('-------------------------------\n');
fprintf('Task 5: With Jacobi iteration method:\n')
task_iteration("Jacobi", 1000000, 1e-7, [8, 16, 32, 64, 128]);

%% TASK 6
fprintf('-------------------------------\n');
fprintf('Task 6: With Gauss-Seidel iteration method:\n');
task_iteration("Gauss-Seidel", 1000000, 1e-7, [8, 16, 32, 64, 128]);

%% TASK 7
fprintf('-------------------------------\n');
fprintf('Task 7: With gradient descent iteration method:\n');
task_iteration("gradient descent", 1000000, 1e-7, [8, 16, 32, 64, 128]);

%% TASK 8
fprintf('-------------------------------\n');
fprintf('Task 8: With conjugate gradient iteration method:\n');
task_iteration("conjugate gradient", 1000000, 1e-7, [8, 16, 32, 64, 128]);

%% local functions
function task_iteration(method, max_iteration, tol, n_array)
    % task_iteration - using iteration method to iterate with n_array and check their convergence.

    method = string(method);
    for n = n_array
        fprintf('For n=%d:\n', n);
        h = 1 / n;
        A = csr_tri_diag_matrix(n);
        b = h * ones(n, 1);
        x0 = zeros(n, 1);
        if method == "conjugate gradient"
            cg_solver(A, b, x0, max_iteration, tol);
        else
            iteration_solver(A, b, x0, max_iteration, tol, method);
        end
    end
end
