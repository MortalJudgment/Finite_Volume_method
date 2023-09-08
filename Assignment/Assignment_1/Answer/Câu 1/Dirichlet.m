
% Finite Volume method for Laplace Equation on 1D
% Consider on Domain [0 1]
% with Dirichlet condition - homogeneous
function u = Dirichlet(x,x_cp,N,u0,uN,v)
%-------------------------------------------------------------------------%
%----- Solve 1D Laplace equation: -----%
%                 -u_xx = f(x)          ,in [a,b]
%-------------------------------------------------------------------------%
% Creare the Matrix
A=zeros(N+2,N+2);
for ii=1:N+2
    if(ii==1)
        A(ii,ii)   = 1;
    elseif(ii==N+2)
        A(ii,ii)   = 1;
    else
        alpha_i     = -1.0/((x(ii)-x(ii-1))*(x_cp(ii)-x_cp(ii-1)));
        gamma_i     = -1.0/((x(ii)-x(ii-1))*(x_cp(ii+1 )-x_cp(ii)));
        A(ii,ii-1) = alpha_i - v/(2*(x(ii)-x(ii-1)));
        A(ii,ii)   = -(alpha_i+gamma_i);
        A(ii,ii+1) = gamma_i + v/(2*(x(ii)-x(ii-1)));
    end
end

% Create vector b
b=zeros(N+2,1);
for ii=1:N+2
    if(ii==1)
        b(ii)   = u0;
    elseif(ii==N+2)
        b(ii)   = uN;
    else
%         b(ii)  =   f(1/2.0*x(ii) + 1/2.0*x(ii-1));       % Midpoint rule
        b(ii)   =   1/2.0*f(x(ii)) + 1/2.0*f(x(ii-1));    % Trepozoidal rule
    end
end

% u=zeros(N,1);
u = A\b;