function u_ex=exact_solution_u(x,y,k);

if (k==1)
    u_ex=(1-cos(2*pi*x))*sin(2*pi*y);
end

if (k==2)
   u_ex=- cos(2*pi*x)*sin(2*pi*y) + sin(2*pi*y); 
end
if (k==3)
    u_ex= x^2*(x-1)^2*(4*y^3-6*y^2+2*y);
end