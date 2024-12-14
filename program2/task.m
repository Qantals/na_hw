clear; clc;
%% Task 1
title_s = "Task 1: 3 intepolation approximate 1st derivation";
fprintf(title_s + "\n")
f = @(x) exp(x);
f_d1 = @(x) exp(x);
x = 0;
h = logspace(-1, -8, 8);

approx = d1_pol(f, x, h);
gtruth = f_d1(x);
format shortE
error = abs(gtruth - approx)

figure;
loglog(h, error, '-o');
title(title_s);
% saveas(gcf, 'task1.png');

%% Task 2
fprintf('-------------------------------\n');
fprintf("Task 2: composite integration with 'midpoint', 'trapezoidal', 'Simpson' and 'Gauss-3points'\n");
a = 0;
b = 1;
n = [1,2,4,8,16,32]
f =@(x) x .* exp(x);
gtruth = f(b) - f(a) - (exp(b) - exp(a));
format shortE

fprintf('********* midpoint method: *********\n');
mid = int_comp(f,a,b,n,'midpoint');
error_mid = abs(mid - gtruth)

fprintf('********* trapezoidal method: *********\n');
trap = int_comp(f,a,b,n,'trapezoidal');
error_trap = abs(trap - gtruth)

fprintf('********* Simpson method: *********\n');
simpson = int_comp(f,a,b,n,'Simpson');
error_simpson = abs(simpson - gtruth)

fprintf('********* Gauss-3points method: *********\n');
gauss3 = int_comp(f,a,b,n,'Gauss-3points');
error_gauss3 = abs(gauss3 - gtruth)

%% Task 3
fprintf('-------------------------------\n');
fprintf("Task 3: numerical differential with'Euler-e', 'trapezoidal-e', 'RK-3' or 'RK-4'\n");
yd = @(x,y) y .* (1-y);
y0 = 0.1;
x = 2;
h = [0.1, 0.05, 0.025, 0.0125]
gtruth = (y0 * exp(x)) / (1 - y0 + y0 * exp(x));
format shortE

fprintf('********* Euler-e method: *********\n');
y_euler = diff(yd,y0,x,h,'Euler-e');
error_euler = abs(y_euler - gtruth)

fprintf('********* trapezoidal-e method: *********\n');
y_trape = diff(yd,y0,x,h,'trapezoidal-e');
error_trape = abs(y_trape - gtruth)

fprintf('********* RK-3 method: *********\n');
y_rk3 = diff(yd,y0,x,h,'RK-3');
error_rk3 = abs(y_rk3 - gtruth)

fprintf('********* RK-4 method: *********\n');
y_rk4 = diff(yd,y0,x,h,'RK-4');
error_rk4 = abs(y_rk4 - gtruth)

