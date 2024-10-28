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
csr_jacobi_iteration(A, b, b)

%% TASK 4
fprintf('Task 4: csr_gs_iteration(A, b, b)\n');
csr_gs_iteration(A, b, b)

%% TASK 5
fprintf('-------------------------------\n');
fprintf('Task 5: With Jacobi iteration method:\n')
max_iteration = 1000000;
tol = 1e-7;
for n = [8, 16, 32, 64, 128]
    fprintf('For n=%d:\n', n);
    h = 1 / n;
    A = csr_tri_diag_matrix(n);
    b = h * ones(n, 1);
    x0 = zeros(n, 1);
    x = jacobi_solver(A, b, x0, max_iteration, tol);
end

%% TASK 6
fprintf('-------------------------------\n');
fprintf('Task 6: With Gauss-Seidel iteration method:\n');
max_iteration = 1000000;
tol = 1e-7;
for n = [8, 16, 32, 64, 128]
    fprintf('For n=%d:\n', n);
    h = 1 / n;
    A = csr_tri_diag_matrix(n);
    b = h * ones(n, 1);
    x0 = zeros(n, 1);
    x = gs_solver(A, b, x0, max_iteration, tol);
end
