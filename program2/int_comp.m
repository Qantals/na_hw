function res = int_comp(f, a, b, n, method)
    % numerical composite intergration with 4 methods.
    %
    % Inputs:
    %     f: function handle
    %     a,b: int, integration interval, like a=0,b=1
    %     n: int | vector, # composite intervals
    %     method: str, Choose 'midpoint', 'trapezoidal', 'Simpson' or 'Gauss-3points'
    %
    % Outputs:
    %     res: float | vector, approximate integration

    switch method
        case "midpoint"
            int_fn = @(f,h,i,c,d) 2*d * f(c);
        case "trapezoidal"
            int_fn = @(f,h,i,c,d) d * (f(h(i)) + f(h(i+1))) ;
        case "Simpson"
            int_fn = @(f,h,i,c,d) d/3 * (f(h(i)) + 4*f(c) + f(h(i+1)));
        case "Gauss-3points"
            int_fn = @(f,h,i,c,d) d * (5/9*f(c - sqrt(3/5)*d) + 8/9*f(c) + 5/9*f(c + sqrt(3/5)*d));
        otherwise
            error("Invalid method. Choose 'midpoint', 'trapezoidal', 'Simpson' or 'Gauss-3points'.");
    end

    res = zeros(1, length(n));
    for label_n = 1:length(n)
        h = linspace(a,b,n(label_n)+1);
        for i = 1:n(label_n)
            c = (h(i+1) + h(i)) / 2;
            d = (h(i+1) - h(i)) / 2;
            res(label_n) = res(label_n) + int_fn(f,h,i,c,d);
        end
    end

end
