function yarr = diff(yd, y0, xend, harr, method)
    % numerical differential with 4 methods.
    %
    % Inputs:
    %     yd(x,y): function handle, first derivation of y.
    %     y0: float, y(0)
    %     xend: float, res = y(xend)
    %     harr: float | vector, step size for x
    %     method: str, Choose 'Euler-e', 'trapezoidal-e', 'RK-3' or 'RK-4' ('-e' means explicity method)
    %
    % Outputs:
    %     yarr = y(xend): float | vector, approximate differential

    yarr = zeros(1, length(harr));
    for label_h = 1:length(harr)
        h = harr(label_h);
        x = 0;
        y = y0;
        while (x < xend)
            switch method
                case "Euler-e"
                    y = y + h * yd(x,y);
                case "trapezoidal-e"
                    y = y + h/2 * (yd(x,y) + yd(x + h, y + h*yd(x,y)));
                case "RK-3"
                    k1 = h * yd(x,y);
                    k2 = h * yd(x + h/2, y + k1/2);
                    k3 = h * yd(x + h, y - k1 + 2*k2);
                    y = y + (k1 + 4*k2 + k3) / 6;
                case "RK-4"
                    k1 = h * yd(x,y);
                    k2 = h * yd(x + h/2, y + k1/2);
                    k3 = h * yd(x + h/2, y + k2/2);
                    k4 = h * yd(x + h, y + k3);
                    y = y + (k1 + 2*k2 + 2*k3 + k4) / 6;
                otherwise
                    error("Invalid method. Choose 'Euler-e', 'trapezoidal-e', 'RK-3' or 'RK-4'.");
            end
            x = x + h;
        end
        yarr(label_h) = y;
    end

end
