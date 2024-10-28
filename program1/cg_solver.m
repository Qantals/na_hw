function x = cg_solver(A, b, x0, max_iteration, tol)
    % cg_solver - Congugate gradient iteration method solver.

    x = x0;
    r = b - csr_vmult(A, x);
    p = r;
    for k = 1 : max_iteration
        if norm(r) < tol || (p' * csr_vmult(A, p)) < tol
            fprintf('Iteration stops at step %d.\n', k);
            return;
        end
        [x, r, p] = csr_cg_iteration(A, x, r, p);
    end
    fprintf('Max iteration step reached.\n');
end
