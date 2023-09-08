%-------------------------------------------------------------------------%
%----- Solve 1D Laplace equation: -----%
%                 -u_xx = f(x)          ,in [a,b]
%-------------------------------------------------------------------------%
clc
clear all
close all
%----------------%
% Domain of x
ax = 0.0;
bx = 1.0;
%Coefficience
beta = 1/16;
%----------------%
N = 50 ;                      % Number of control volume
%------------------------------------------------------------------------%

dx = (bx-ax)/N;
dt = 1/2*1/(2*beta)*dx^2;
% Create the mesh point
x = zeros(N+1,1);
for ii=1:N+1
    x(ii) = ax+(ii-1)*dx;                 % uniform_Mesh
    %     x(ii) = bx-cos(pi/2*(i_iter-1)/N);    non-uniform_Mesh
end

% Create control point
x_cp = zeros(N+2,1);
for ii=1:N+2
    if(ii==1)
        x_cp(ii) = x(ii);
    else
        if(ii==N+2)
            x_cp(ii) = x(ii-1);
        else
            x_cp(ii) = (1*x(ii-1)+1*x(ii))/2.0;
        end
    end
end
%-----------------%
%Initial condition
u1 = (u0(x_cp));
%------------------------------------------------------------------------%
for jj=1:600
    clf    
    t_k = jj*dt;
%---------- Dicrete solution -----------%
    u_dis = Dirichlet(beta,N,dt,x,x_cp,u1);
    u1 = u_dis;
%----------- Exact solution ------------%
    u_ex = zeros(N+2,1);
    for ii=1:N+2
        u_ex(ii) = uex(x_cp(ii),t_k);
    end
%---------------------------------------%
%-------------- Drawing ----------------%
%---------------------------------------%
    hold on
    plot(x_cp,u_dis,'or',x_cp,u_ex,'b');
    axis([0 1 -1 1])
    pause(10^(-6))
end
hold off