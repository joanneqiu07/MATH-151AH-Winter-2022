% the three definite integrals in problem 4
f = @(x) x.*log(x);
g = @(x) x./(x.^2 + 4);
h = @(x) 1/sqrt(2* pi) * exp(-(x.^2)/2);

% the fourth Legendre polynomial
L = @(x) x^4 - (6/7)*x^2 + 3/35;

% intialize a matrix
x = zeros(4,1);

% use build-in fucntion fzero() to find the roots
x(1) = fzero(L,-0.9);
x(2) = fzero(L,-0.3);
x(3) = fzero(L,0.3);
x(4) = fzero(L,0.9);

% calculate the weights
w = weights(x,4);

% evaluate the definite integral
disp("The integral of 4a:")
gauss4(f,x,1,2)

function p = P(n,t)

if n == 0
    p = 1;
    return
end

if n == 1
    p = t;
    return
end

temp1 = 1;
temp2 = t;
temp3 = 0;

for i = 1:(n-1)
    temp3 = temp2;
    temp2 = ((2*i + 1)*t*temp2 - i*temp1)/(i+1);
    temp1 = temp3;
end

p = temp2;
end

function w = weights(x,n)

Pol = zeros(n);
b = zeros(n,1);
b(1,1) = 2;

for i = 1:n
    for j = 1:n
        if n==1
            Pol(i,j) = P(i-1,x);
        else
            Pol(i,j) = P(i-1,x(j));
        end
    end
end

w = Pol\b;
end

function gauss = gauss4(f,x,a,b)

xs = (b-a)/2*x + (a+b)/2;
vals = f(xs);

w = weights(x,4);

gauss = (b-a)/2 * dot(w,vals);

end