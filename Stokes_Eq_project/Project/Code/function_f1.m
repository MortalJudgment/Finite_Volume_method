function f_1=function_f1(x,y,k);

if (k==1)
    f_1=y + 3*x^2*y^2 - 4*pi^2*cos(2*pi*x)*sin(2*pi*y) - 4*pi^2*sin(2*pi*y)*(cos(2*pi*x) - 1) + 1;
end

if (k==2)
    f_1=4*pi^2*sin(2*pi*x) + 4*pi^2*sin(2*pi*y) - 8*pi^2*cos(2*pi*x)*sin(2*pi*y);
end
if (k==3)
    f_1= 2*x - 2*(x-1)^2*(4*y^3-6*y^2+2*y) - 2*x^2*(4*y^3-6*y^2+2*y) - x^2*(24*y-12)*(x-1)^2 - 4*x*(2*x-2)*(4*y^3-6*y^2+2*y);

end
