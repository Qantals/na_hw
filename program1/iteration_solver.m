function x = iteration_solver(A, b, x0, max_iteration, tol, method)
    % iteration_solver - Iteration method solver.

    x = x0;
    for k = 1 : max_iteration
        p = csr_iteration(A, b, x, method);
        if norm(x - p) < tol
            fprintf('Iteration stops at step %d.\n', k);
            x = p;
            return;
        end
        x = p;
    end
    fprintf('Max iteration step reached.\n');
end
