function r = uex(x,y,test)
if nargin == 2
    test = 1;
end
switch test
    case 1
        r(1) = (1 - cos(2*pi*x))*sin(2*pi*y);
        r(2) = -(1 - cos(2*pi*x))*sin(2*pi*y);
    case 2
        r(1) = (1 - cos(2*pi*x))*sin(2*pi*y);
        r(2) = -(1 - cos(2*pi*y))*sin(2*pi*x);
    case 3
        r(1) = x^2*(x-1)^2*(4*y^3-6*y^2+2*y);
        r(2) = - y^2*(y-1)^2*(4*x^3-6*x^2+2*x);
    otherwise
        fprintf('Undentified valuable!');
end