function r = pex(x,y,test)
if nargin == 2
    test = 1;
end
switch test
    case 1
        r = x*y+x+y+x^3*y^2 - 4.0/3.0;
    case 2
        r = -2*pi*(cos(2*pi*x) - cos(2*pi*y));
    case 3
        r = x^2 + y^2 - 2.0/3.0;
    otherwise
        fprintf('Undentified valuable!');
end