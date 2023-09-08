% Finite Volume method for Laplace Equation on 1D
% Consider on Domain [0 1]
% with Dirichlet condition - homogeneous
function u = Neumann(x,x_cp,dx,N,du0,duN)
%-------------------------------------------------------------------------%
%----- Solve 1D Laplace equation: -----%
%                 -u_xx = f(x)          ,in [a,b]
%-------------------------------------------------------------------------%
% Creare the Matrix
A=zeros(N+2,N+2);
for ii=1:N+2
    if(ii==1)
        A(ii,ii)   = -1/(x_cp(ii+1)-x_cp(ii));
        A(ii,ii+1)   = 1/(x_cp(ii+1)-x_cp(ii));
    elseif(ii==N+2)
        for iz = 1:N
            A(ii,iz+1) = x(iz+1)-x(iz);
        end
    else
        alpha_i     = -1.0/((x(ii)-x(ii-1))*(x_cp(ii)-x_cp(ii-1)));
        gamma_i     = -1.0/((x(ii)-x(ii-1))*(x_cp(ii+1)-x_cp(ii)));
        A(ii,ii-1) = alpha_i;
        A(ii,ii)   = -(alpha_i+gamma_i);
        A(ii,ii+1) = gamma_i;
    end
end

% Create vector b
b=zeros(N+2,1);
for ii=1:N+2
    if(ii==1)
        b(ii)   = du0;
    elseif(ii==N+2)
        b(ii)   = duN;
    else
    %        b(ii)  =   f((x(ii) + x(ii+1))/2.0);       % Midpoint rule
        b(ii)   =   (f(x(ii)) + f(x(ii-1)))/2.0;    % Trepozoidal rule
    end
end
epsilon11=0.0000001;
% u=zeros(N,1);
u = (A+epsilon11*eye(N+2,N+2))\b;