% Test the exp function.
fun = inline('sin(x)');
fun_der = inline('cos(x)');
interval = [-8, 0];
order = 5;

A = remez(fun, fun_der, interval, order);

A1 = A(1:end-1);  % the 3 coefficients of the second order polynomial
E = A(end); % the maximum approximation error 

numTestPoints = 100;
x = linspace(interval(1), interval(2), numTestPoints);
fx = feval(fun, x);
fest = polyval(A1(end:-1:1), x-interval(1));
ferr = err(x, fun, A1, interval(1));
