a = input("a = ");
b = input("b = ");
m = input("m = ");

fun = input('f(x) = ', 's');
f = str2func(['@(x) ' fun ]);

out = simps(a,b,m,f);

function out = simps(a, b, m, f)
x = zeros();
out = 0;
h = (b - a)/(2*m);
    for i = 0:2*m
        x(i+1) = a + i*h;
    end

    for j = 1:m
        out = out + h/3*(f(x(2*j-1))+4*f(x(2*j))+f(x(2*j+1)));
    end
end