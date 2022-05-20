% the matrix A
A = [2.4117, 0.6557, 0.6787, 0.6555; 0.9157, 1.8804, 0.7577, 0.1712; 
    0.7922, 0.8491, 3.0905, 0.7060; 0.9595, 0.9340, 0.3922, 2.3175];
m = 4;

% the matrix b
b = [8.3813; 7.6345; 14.5862; 13.2743];

% other variables
x = rand(m, 1);
y = A\b;
err = norm(x - y)/norm(y);
tol = 10^(-12);

% Jacobi's method
D = diag(diag(A));
Di = inv(D);
J = A - D;

% Gauss-Seidel method
N = [2.4117, 0, 0, 0; 0.9157, 1.8804, 0, 0; 
    0.7922, 0.8491, 3.0905, 0; 0.9595, 0.9340, 0.3922, 2.3175];
Ni = inv(N);
P = A - N;

% SOR method
w = 0.5;
A2 = D + w*(N - D);
A2i = inv(A2);
A1 = (w - 1)*D + w*P;

[s,c,i] = sor(x, y, err, tol, A2i, A1, w, b);
[j,a,k] = jac(x, y, err, tol, Di, J, b);
[g,d,l] = gau(x, y, err, tol, Ni, P, b);

iter1 = 1:i;
iter2 = 1:k;
iter3 = 1:l;

plot(iter1,log(c(iter1)),'Color','b','Marker','o')
hold on
plot(iter2,log(a(iter2)),'Color','m','Marker','o')
hold on
plot(iter3,log(d(iter3)),'Color','g','Marker','o')

% function of Jacobi's method
function [j,a,k] = jac(x, y, err, tol, Di, J, b)
j = x;
a = zeros();
a(1) = err;
k = 1;
    while err > tol
        k = k + 1;
        j = Di * (b - J*j);
        err = norm(j - y)/norm(y);
        a(k) = err;
    end
end

% function of Gauss-Seidel's method
function [g,d,l] = gau(x, y, err, tol, Ni, P, b)
g = x;
d = zeros();
d(1) = err;
l = 1;
    while err > tol
       l = l+1;
       g = Ni * (b - P*g);
       err = norm(g - y)/norm(y);
       d(l) = err;
    end
end

% function of SOR method
function [s,c,i] = sor(x, y, err, tol, A2i, A1, w, b)
s = x;
c = zeros();
c(1) = err;
i = 1;
    while err > tol
        i = i + 1;
        s = A2i*(-A1*s) + A2i*w*b;
        err = norm(s - y)/norm(y);
        c(i) = err;
    end
end
