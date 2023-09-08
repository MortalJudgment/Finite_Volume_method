function r = g(x,y,test)
if nargin == 2
    test = 1;
end
switch test
    case 1
        r = x + 2*x^3*y + 4*pi^2*cos(2*pi*y)*sin(2*pi*x) + 4*pi^2*sin(2*pi*x)*(cos(2*pi*y) - 1) + 1;
    case 2
        r = 4*pi^2*cos(2*pi*y)*sin(2*pi*x) - 4*pi^2*sin(2*pi*y) + 4*pi^2*sin(2*pi*x)*(cos(2*pi*y) - 1);
    case 3
        r = 2*y + 2*(y-1)^2*(4*x^3-6*x^2+2*x) + 2*y^2*(4*x^3-6*x^2+2*x) + y^2*(24*x-12)*(y-1)^2 + 4*y*(2*y-2)*(4*x^3-6*x^2+2*x);
    otherwise
        fprintf('Undentified valuable!');
end