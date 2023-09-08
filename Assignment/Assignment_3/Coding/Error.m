% Finite Volume method for Diffusion Advection Equation on 1D
%-------------------------------------------------------------------------%
%----- Solve 1D Diffusion Advection Equation: ----------------------------%
%----         u_t + alpha*u_x = beta*u_xx      ,in Omega              ----%
%-------------------------------------------------------------------------%
% Consider on Domain Omega = (0 1), t>0
% with initial condition u0(x)
% and Dirichlet boudary condition u(0,t)=g1(t); u(1,t)=g2(t)
clc
close all
clear all
format long
%-------------------------------------------------------------------------%
%----------------%
% Domain of x
ax = 0.0;
bx = 1.0;
%Coefficience
alpha = 0;          %Coefficience Advection term
beta = 1/16;        %Coefficience Diffusion term
%----------------%
Nx = 20 ;                      % Number of control volume
%------------------------------------------------------------------------%
m = 4;
h = zeros(m,1);
errorL2 = zeros(m,1);
error2 = zeros(m,1);
error3 = zeros(m,1);
%------------------------------------------%
T = 0.1; % number of step on time
for q = 1:m
    h(q) = 1/Nx;
    k = 1/2*1/(2*beta/h(q)^2+alpha/h(q));
    iter = round(T/k);
    % Create the mesh point
    x = zeros(Nx+1,1);
    for ii=1:Nx+1
        x(ii) = ax+(ii-1)*h(q);                           % uniform_Mesh
    end
    
    % Create control point
    x_cp = zeros(Nx+2,1);
    for ii=1:Nx+2
        if(ii==1)
            x_cp(ii) = x(ii);
        else
            if(ii==Nx+2)
                x_cp(ii) = x(ii-1);
            else
                x_cp(ii) = (1*x(ii-1)+1*x(ii))/2.0;
            end
        end
    end
    u1 = u0(x_cp);
    u2 = u0(x_cp);
    u3 = u0(x_cp);
    Temp = 0;
    for jj=1:iter
        t_k = jj*k;
%         jj
        u1 = Dirichlet(alpha,beta,1,Nx,k,x,x_cp,u1);        %Backward Euler
%         u2 = Dirichlet(alpha,beta,1/2,Nx,k,x,x_cp,u2);      %Crank-Nilcoson
%         u3 = Dirichlet(alpha,beta,0,Nx,k,x,x_cp,u3);        %Forward Euler
    end
    e_ii = zeros(Nx+2,1);
    for ii=1:Nx+2
        e_ii(ii) = abs(u1(ii)-uex(t_k,x_cp(ii),alpha));
        Temp = Temp + e_ii(ii)^2*h(q);
        %             temp2(jj) = temp2(jj)+abs(u2(i)-uex(t_k,x_cp(i),alpha))^2*h(q);
        %             temp3(jj) = temp3(jj)+abs(u3(i)-uex(t_k,x_cp(i),alpha))^2*h(q);
    end
    errorL2(q) = sqrt(Temp);
%     error2(q) = sqrt((temp2(iter)));
%     error3(q) = sqrt((temp3(iter)));
    Nx = Nx*2;
end
plot(-log(h),-2*log(h)-3,-log(h),-1*log(h)-1,-log(h),-log(errorL2)+2)
legend('Oder 2','Oder 1','Error Backward Euler Method')
title('Order of Error (consider spaceline) for each method')
%------------------------------------------%