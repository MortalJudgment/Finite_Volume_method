% Finite Volume method for Laplace Equation on 1D
% Consider on Domain [0 1]
% with Dirichlet condition - homogeneous
function u = Dirichlet(beta,N,dt,x,x_cp,u1)
%-------------------------------------------------------------------------%
%----- Solve 1D Laplace equation: -----%
%                 u_t = beta*u_xx          ,in [a,b]
%-------------------------------------------------------------------------%

% Creare the Matrix
A=zeros(N+2,N+2);
for ii=1:N+2
    if(ii==1)
        A(ii,ii)    = 1;
    elseif(ii==N+2)
        A(ii,ii)    = 1;
    else
        alpha_i     = -1.0/((x(ii)-x(ii-1))*(x_cp(ii)-x_cp(ii-1)));
        gamma_i     = -1.0/((x(ii)-x(ii-1))*(x_cp(ii+1)-x_cp(ii)));
        A(ii,ii-1)  = -alpha_i;
        A(ii,ii)    = (alpha_i+gamma_i);
        A(ii,ii+1)  = -gamma_i;
    end
end

% Create vector b
b=zeros(N+2,1);
for ii=1:N
    if(ii==1)
        b(ii)   = 0;
    elseif(ii==N+2)
        b(ii)   = 0;
    else
    %        b(ii)  =   f((x(ii) + x(ii+1))/2.0);       % Midpoint rule
        b(ii)   =   (f(x(ii)) + f(x(ii-1)))/2.0;    % Trepozoidal rule
    end
end

u = u1 + dt*beta*(A*u1+b);