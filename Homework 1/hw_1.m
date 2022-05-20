% Problem 9
% Declare function handle
% f = @(x) 10*x.^2 - 4*x - 5;
% df = @(x) 20*x - 4;
% f = @(x) x.^4 - x.^3 - 3*x.^2 + 5*x - 2;
% df = @(x) 4*x.^3 - 3*x.^2 - 6*x + 5;
% f = @(x) sin(x) - log(x);
% df = @(x) cos(x) - 1/x;
g = @(x) x.^2 + 4*cos(x);
f = @(x) 2*x - 4*sin(x);
df = @(x) 2 - 4*cos(x);

% initial values
x_0 = 1; 
% x_0 = 1.3;
% x_0 = 2;

% choose x_1 for the secant method
% x_1 = x_0 - 0.1;
% x_1 = x_0 + 0.1;

% alpha value
% alpha = -0.1;
% alpha = 0.1;

% Max number of iteration
max_iter = 4;


% [fix,x] = fixer(f,x_0,alpha,max_iter);
% f(fix)
% disp(fix)
[newt,y] = newton(f,df,x_0,max_iter);
f(newt)
disp(newt)
% [sec,z] = secant(f,x_0,x_1,max_iter);
% [steff,w] = steffen(f,x_0,max_iter);
% plot(x,f(x),'Color','b','Marker','o')
% hold on
fplot(g, [1 2], 'Color', 'b', 'LineStyle','-');
% plot(y,f(y),'Color','r','Marker','o')
% hold on
% plot(z,f(z),'Color','g','Marker','o')
% hold on
% plot(w,f(w),'Color','m','Marker','o')


function [fix,x] = fixer(f,x_0,alpha,max_iter)
fix = x_0;
x = zeros(1,max_iter);
% Run a for loop
    for i = 1:max_iter
        fix = fix - alpha*f(fix);
        x(i) = fix;
    end
end

function [newt,y] = newton(f,df,x_0,max_iter)
newt = x_0;
y = zeros(1,max_iter);
% Run a for loop
    for i = 1:max_iter
        newt = newt -f(newt)/df(newt);
        y(i) = newt;
    end
end

function [sec,z] = secant(f,x_0,x_1,max_iter)
sec = x_1;
sec_0 = x_0;
z = zeros(1,max_iter+1);
z(1) = x_1;
% Run a for loop
    for i = 1:max_iter
        sec = sec - (f(sec)*(sec - sec_0))/(f(sec)-f(sec_0));
        z(i+1) = sec;
        sec_0 = z(i);
    end
end

function [stef,w] = steffen(f,x_0,max_iter)
stef = x_0;
w = zeros(1,max_iter);
% Run a for loop
    for i = 1:max_iter
        stef = stef - (f(stef)^2)/(f(stef+f(stef))-f(stef));
        w(i) = stef;
    end
end

