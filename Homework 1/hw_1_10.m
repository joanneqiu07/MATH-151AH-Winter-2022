% Problem 10
% declare function
% df of 10*x.^2 - 4*x + 5
f = @(x) 20*x - 4;
% df of f is 20

% the tolerance
tol = 0.00001;

% chose an initial value: -1.5, -1, -0.7, -0.2, -0.1, 0.3, 0.5, 0.8, 1, 1.4
x_0 = 1.4;


% [fix1,x1] = fixer(f,x_0,0.01,tol);
% [fix2,x2] = fixer(f,x_0,0.05,tol);
% [fix3,x3] = fixer(f,x_0,2,tol);
% uncomment the following lines to  display the results
% disp("x1 w/ alpha = 0.01:");
% fprintf('%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n', x1);
% fprintf('\n\n')
% disp("x2 w/ alpha = 0.05:");
% fprintf('%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f', x2);
% fprintf('\n\n')
% disp("x3 w/ alpha = 2:");
% fprintf('%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n', x3);
% fprintf('\n')

[newt,y] = newton(f,x_0,tol);
% uncomment the following lines to display the results
disp("newton's method finding root:");
fprintf('%9.6f %9.6f %9.6f %9.6f %9.6f\n', y);
fprintf('\n')

function [fix,x] = fixer(f,x_0,alpha,tol)
fix = x_0;
x = zeros();
x(1) = x_0;
fix = fix - alpha*f(fix);
x(2) = fix;
i = 2;
% Run a while loop
    while(abs(x(i)- x(i-1))>=tol)
        i = i + 1;
        fix = fix - alpha*f(fix);
        x(i) = fix;
    end
end

function [newt,y] = newton(f,x_0,tol)
newt = x_0;
y = zeros();
y(1) = x_0;
newt = newt - f(newt)/20;
y(2) = newt;
j = 2;
% Run a while loop
    while(abs(y(j)- y(j-1)) >= tol)
        j = j + 1;
        newt = newt -f(newt)/20;
        y(j) = newt;
    end
end