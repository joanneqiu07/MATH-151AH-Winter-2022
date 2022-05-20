% problem 3
f1 = @(x,y) x.^2 + y.^2 + 5*x;
f2 = @(x,y) 2*x*y + 3*y.^2 + y;
j1 = @(x) 2*x + 5;
j2 = @(y) 2*y;
j3 = @(y) 2*y;
j4 = @(x,y) 2*x + 6*y + 1;

% x^(0)
x_0 = -1;
y_0 = -1;

% the number of iteration for the Jacobi's method
iter1 = 200;
iter2 = 100;
iter3 = 80;

[newt1, newt2,x,y,z] = newton(f1,f2,x_0,y_0,j,j1,j2,j3,j4,iter1);
[newt3, newt4,m,n,l] = newton(f1,f2,x_0,y_0,j,j1,j2,j3,j4,iter2);
[newt5, newt6,a,b,c] = newton(f1,f2,x_0,y_0,j,j1,j2,j3,j4,iter3);
plot(x,y,'Color','b','Marker','o')
hold on
plot(m,n,'Color','c','Marker','o')
hold on
plot(a,b,'Color','m','Marker','o')
hold on

function j = jac(x, A, b, max_iter)
D = diag(diag(A));
Di = inv(D);
J = A - D;
j = x;
    for i = 1:max_iter
        j = Di * (b - J*j);
    end
end

function [newt1, newt2,x,y,z] = newton(f1,f2,x_0,y_0,j,j1,j2,j3,j4,max_iter)
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
M = M - jac([x_0; y_0], J, F, max_iter);
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
        M = M - jac([x_0; y_0], J, F, max_iter);
        newt1 = M(1);
        newt2 = M(2);
        x(j) = newt1;
        y(j) = newt2;
        disp(y(j));
        z(j) = norm([x(j);y(j)]);
    end
end