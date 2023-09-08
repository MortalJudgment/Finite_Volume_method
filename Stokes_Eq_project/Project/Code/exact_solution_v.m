function v_ex=exact_solution_v(x,y,k);

if (k==1)
    v_ex=-(1-cos(2*pi*y))*sin(2*pi*x);
end

if(k==2)
    v_ex=cos(2*pi*y)*sin(2*pi*x) - sin(2*pi*x);
end
if(k==3)
    v_ex=- y^2*(y-1)^2*(4*x^3-6*x^2+2*x);
end
