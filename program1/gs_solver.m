function x = gs_solver(A, b, x0, max_iteration, tol)
    x = x0;
    for k = 1 : max_iteration
        p = csr_gs_iteration(A, b, x);
        if norm(x - p) < tol
            fprintf('Iteration stops at step %d.\n', k);
            x = p;
            return;
        end
        x = p;
    end
    fprintf('Max iteration step reached.\n');
end
