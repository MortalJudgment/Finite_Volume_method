% Finite Volume method for Diffusion Advection Equation on 1D
% Consider on Domain Omega = (0 1), t>0
% with initial condition u0(x)
% and Dirichlet boudary condition u(0,t)=g1(t); u(1,t)=g2(t)
function u = Dirichlet(alpha,beta,theta,N,dt,x,x_cp,u1)
%-------------------------------------------------------------------------%
%-----  Solve 1D Laplace equation:  -----%
%                 u_t + alpha*u_x = beta*u_xx   ,in [a,b]
%-------------------------------------------------------------------------%
% unit matrix
I = eye(N+2,N+2);
% Creare the Matrix
A=zeros(N+2,N+2);
for ii=2:N+1
    alpha_i     = 1.0/((x(ii)-x(ii-1))*(x_cp(ii)-x_cp(ii-1)));
    gamma_i     = 1.0/((x(ii)-x(ii-1))*(x_cp(ii+1)-x_cp(ii)));
    A(ii,ii-1)  = beta*dt*theta*alpha_i + alpha*dt/(2*(x(ii)-x(ii-1)));
    A(ii,ii)    = beta*dt*theta*(-alpha_i-gamma_i);
    A(ii,ii+1)  = beta*dt*theta*gamma_i - alpha*dt/(2*(x(ii)-x(ii-1)));
end

% Create vector B
B=zeros(N+2,N+2);
for ii=2:N+1
    alpha_i     = 1.0/((x(ii)-x(ii-1))*(x_cp(ii)-x_cp(ii-1)));
    gamma_i     = 1.0/((x(ii)-x(ii-1))*(x_cp(ii+1)-x_cp(ii)));
    B(ii,ii-1)  = beta*dt*(1-theta)*alpha_i;
    B(ii,ii)    = beta*dt*(1-theta)*(-alpha_i-gamma_i);
    B(ii,ii+1)  = beta*dt*(1-theta)*gamma_i;
end

u = (I-B)^(-1)*(I+A)*u1;