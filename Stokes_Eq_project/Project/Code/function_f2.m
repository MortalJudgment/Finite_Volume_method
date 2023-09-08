function f_2=function_f2(x,y,k);

if (k==1)
    f_2=x + 2*x^3*y + 4*pi^2*cos(2*pi*y)*sin(2*pi*x) + 4*pi^2*sin(2*pi*x)*(cos(2*pi*y) - 1) + 1;
end

if(k==2)
    f_2=8*pi^2*cos(2*pi*y)*sin(2*pi*x) - 4*pi^2*sin(2*pi*y) - 4*pi^2*sin(2*pi*x);
end
if(k==3)
    f_2= 2*y + 2*(y-1)^2*(4*x^3-6*x^2+2*x) + 2*y^2*(4*x^3-6*x^2+2*x) + y^2*(24*x-12)*(y-1)^2 + 4*y*(2*y-2)*(4*x^3-6*x^2+2*x);

end


