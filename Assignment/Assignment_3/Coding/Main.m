% Finite Volume method for Diffusion Advection Equation on 1D
% Consider on Domain Omega = (0 1), t>0
% with initial condition u0(x)
% and Dirichlet boudary condition u(0,t)=g1(t); u(1,t)=g2(t)
%-------------------------------------------------------------------------%
%----- Solve 1D Diffusion Advection Equation: ----------------------------%
%----         u_t + alpha*u_x = beta*u_xx      ,in [a,b]              ----%
%-------------------------------------------------------------------------%
clc
clear all
close all
%----------------%
% Domain of x
ax = 0.0;
bx = 1.0;
%Coefficience
alpha = 2;
beta = 1/16;
%----------------%
N = 50 ;                      % Number of control volume
%------------------------------------------------------------------------%

dx = (bx-ax)/N;
dt = 1/2*1/(2*beta/dx^2+alpha/dx);
% Create the mesh point
x = zeros(N+1,1);
for ii=1:N+1
    x(ii) = ax+(ii-1)*dx;                           % uniform_Mesh
    %     x(ii) = bx-cos(pi/2*(i_iter-1)/N);        % non-uniform_Mesh
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
theta = 1;
%------------------------------------------------------------------------%
figh = figure;

for jj=1:600
    clf    
    t_k = jj*dt;
%---------- Dicrete solution -----------%
    u_dis = Dirichlet(alpha,beta,theta,N,dt,x,x_cp,u1);
    u_dis(1) = g(t_k);
    u_dis(N+2) = g(t_k);
%     u_dis(1) = uex(x_cp(1),t_k,alpha);
%     u_dis(N+2) = uex(x_cp(N+2),t_k,alpha);
    u1 = u_dis;
%----------- Exact solution ------------%
    u_ex = zeros(N+2,1);
    for ii=1:N+2
        u_ex(ii) = uex(x_cp(ii),t_k,alpha);
    end
%---------------------------------------%
%-------------- Drawing ----------------%
%---------------------------------------%
    hold on
    plot(x_cp,u_dis,'or',x_cp,u_ex,'b');
    axis([0 1 -1 1])
    title(['Particle at t = ', num2str(jj*dt),' seconds'])
%     pause(10^(-6))
    movieVector(jj) = getframe(figh, [10 10 520 400]);
end
myWriter = VideoWriter('Heat Equation');
myWriter.FrameRate = 20;
open(myWriter);
writeVideo(myWriter,movieVector);
close(myWriter);
hold off
