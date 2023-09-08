function p_ex=exact_pressure(x,y,k);

if (k==1)
    p_ex=x*y+x+y+x.^3*y.^2-4/3;
end

if (k==2)
   p_ex=-2*pi*(cos(2*pi*x) - cos(2*pi*y)); 
end
if (k==3)
    p_ex=(x^2+y^2-2/3);
end
