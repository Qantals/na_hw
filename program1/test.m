clear all;
format long;

n = 64;
A = csr_tri_diag_matrix(n);
b = ones(n, 1) / n;

max_iteration = 1000000;
tol = 1e-7;
x0 = zeros(n, 1);

x = iteration_solver(A, b, x0, max_iteration, tol, "Jacobi");

plot(linspace(0, 1, n + 2), [0, x', 0]);
