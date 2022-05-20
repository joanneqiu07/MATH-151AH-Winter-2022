


choice = 0;
choice = input("Enter 1 for Newton, 0 for Lagrange: ");

xinput = input('Enter x points in the form of [x1, x2,...]:');
x = strrep(xinput, ',', ' ');
yinput = input("Enter y points in the form of [y1, y2,...]:");
y = strrep(yinput, ',', ' ');

xt = x(1):0.1:x(length(x));

if choice == 1
    a = divdiff(x,y);
    yt = evalnewt(a,x,xt);
elseif choice == 0
    m = input("Your choice of x for the Lagrange interpolation:");
    [yt,L] = lag(m,x,xt,y);
end

figure(1);
plot(xt,yt,'-',x,y,'*');
if choice == 0
    disp(yt);
    disp(L);
end

function a = divdiff(x,y) 
n = length(x);
F = zeros(n,n);
F(:, 1) = y;
for i=1:(n-1)
    for j=1:i
        F(i+1,j+1) = (F(i+1,j)-F(i,j))/(x(i+1)-x(i-j+1));
    end
end

a = diag(F);

end

function yt = evalnewt(a,x,xt)
n = length(x) - 1;
yt = a(n)*ones(size(xt));
for i = n:(-1):1
    yt = a(i)+yt.*(xt - x(i));
end
end

function [yt,L] = lag(m,x,xt,y)
yt = zeros(size(xt));
n = length(x) - 1;
sum = 0;
for i = 1:n
    u = 1;
    l = 1;
    for j = 1:n
        if j ~= i
            u = u * (m - x(j));
            l = l * (x(i) - x(j));
        end
    end
    sum = sum + u / l * y(i);
    yt(i) = sum;
end
L = yt(x == m);
end