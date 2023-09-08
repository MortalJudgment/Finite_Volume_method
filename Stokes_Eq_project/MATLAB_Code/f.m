function r = f(x,y,test)
if nargin == 2
    test = 1;
end
switch test
    case 1
        r = y + 3*x^2*y^2 - 4*pi^2*cos(2*pi*x)*sin(2*pi*y) - 4*pi^2*sin(2*pi*y)*(cos(2*pi*x) - 1) + 1;
    case 2
        r = 4*pi^2*sin(2*pi*x) - 4*pi^2*cos(2*pi*x)*sin(2*pi*y) - 4*pi^2*sin(2*pi*y)*(cos(2*pi*x) - 1);
    case 3
        r = 2*x - 2*(x-1)^2*(4*y^3-6*y^2+2*y) - 2*x^2*(4*y^3-6*y^2+2*y) - x^2*(24*y-12)*(x-1)^2 - 4*x*(2*x-2)*(4*y^3-6*y^2+2*y);
    otherwise
        fprintf('Undentified valuable!');
end