function res = d1_pol(f, x, h)
    % Use a 3 points intepolation method to approximate first derivation.
    % Formula: f'(x)=\frac{f(x-2h)-4f(x-h)+3f(x)}{2h}
    %
    % Inputs:
    %     f: function handle
    %     x: float scalar, input variable
    %     h: float scalar | vector, step size
    %
    % Outputs:
    %      res: float scalar | vector, approximation f'(x)

    res = (f(x-2*h)-4*f(x-h)+3*f(x)) ./ (2*h);

end
