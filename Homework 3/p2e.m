% problem 2e
% declare function handle
f1 = @(x,y) x.^2 + y.^2 + 5*x;
f2 = @(x,y) 2*x*y + 3*y.^2 + y;
j1 = @(x) 2*x + 5;
j2 = @(y) 2*y;
j3 = @(y) 2*y;
j4 = @(x,y) 2*x + 6*y + 1;

% tolerance
tol = 0.00001;

% alpha, x_0, y_0
alpha = -0.1;
x_0 = -1;
y_0 = -1;

[fix1, fix2, x, y,z] = fixer(f1, f2,x_0,y_0,alpha,tol);
[newt1, newt2,m,n,l] = newton(f1,f2,x_0,y_0,j1,j2,j3,j4);

iter1 = 1:length(x);
iter2 = 1:50;

%plot(iter1,z(iter1),'Color','b','Marker','o')
%hold on
%plot(iter2,l(iter2),'Color','m','Marker','o')
%hold on
plot(x,y,'Color','b','Marker','o')
hold on
plot(m,n,'Color','m','Marker','o')
hold on

function [fix1, fix2, x, y, z] = fixer(f1, f2,x_0,y_0,alpha,tol)
fix1 = x_0;
fix2 = y_0;
x = zeros();
y = zeros();
z = zeros();
x(1) = x_0;
y(1) = y_0;
z(1) = norm([x(1);y(1)]);
fix1 = fix1 - alpha*f1(fix1,fix2);
fix2 = fix2 - alpha*f2(fix1,fix2);
x(2) = fix1;
y(2) = fix2;
z(2) = norm([x(2);y(2)]);
i = 2;
% Run a while loop
    while(abs(x(i)- x(i-1))>=tol || abs(y(i)- y(i-1))>=tol)
        i = i + 1;
        fix1 = fix1 - alpha*f1(fix1,fix2);
        fix2 = fix2 - alpha*f2(fix1,fix2);
        x(i) = fix1;
        y(i) = fix2;
        z(i) = norm([x(i);y(i)]);
    end
end

function [newt1, newt2,x,y,z] = newton(f1,f2,x_0,y_0,j1,j2,j3,j4)
newt1 = x_0;
newt2 = y_0;
x = zeros();
y = zeros();
z = zeros();
x(1) = x_0;
y(1) = y_0;
z(1) = norm([x(1);y(1)]);
J = [j1(x_0),j2(y_0);j3(y_0),j4(x_0,y_0)];
F = [f1(x_0,y_0);f2(x_0,y_0)];
M = [newt1, newt2];
M = M - J\F;
newt1 = M(1);
newt2 = M(2);
x(2) = newt1;
y(2) = newt2;
z(2) = norm([x(2);y(2)]);
j = 2;
% Run a while loop
    while(abs(x(j)- x(j-1)) >= 0.01 || abs(y(j)- y(j-1))>=0.01)
        j = j + 1;
        J = [j1(newt1),j2(newt2);j3(newt2),j4(newt1,newt2)];
        F = [f1(newt1,newt2);f2(newt1,newt2)];
        M = [newt1, newt2];
        M = M - J\F;
        newt1 = M(1);
        newt2 = M(2);
        x(j) = newt1;
        y(j) = newt2;
        disp(y(j));
        z(j) = norm([x(j);y(j)]);
    end
end