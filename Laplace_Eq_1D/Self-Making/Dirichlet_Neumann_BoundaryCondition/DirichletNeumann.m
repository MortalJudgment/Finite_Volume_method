% Finite Volume method for Laplace Equation on 1D
% Consider on Domain [0 1]
% with Dirichlet condition - homogeneous
function u = DirichletNeumann(ax,bx,N,u0,uN)
%-------------------------------------------------------------------------%
%----- Solve 1D Laplace equation: -----%
%                 -u_xx = f(x)          ,in [a,b]
%-------------------------------------------------------------------------%

dx = (bx-ax)/N;

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
    elseif(ii==N+2)
        x_cp(ii) = x(ii-1);
    else
        x_cp(ii) = (1*x(ii-1)+1*x(ii))/2.0;
    end
end

% Creare the Matrix
A=zeros(N+2,N+2);
for ii=1:N+2
    if(ii==1)
        A(ii,ii)   = 1;
    elseif(ii==N+2)
        A(ii,ii-1) = -2.0/((x(ii)-x(ii-1))*(x_cp(ii)-x_cp(ii-1)));
        A(ii,ii)   = 2.0/((x(ii)-x(ii-1))*(x_cp(ii)-x_cp(ii-1)));
    else
        alpha_i    = -1.0/((x(ii+1)-x(ii))*(x_cp(ii+1)-x_cp(ii)));
        gamma_i    = -1.0/((x(ii+1)-x(ii))*(x_cp(ii+1)-x_cp(ii)));
        A(ii,ii-1) = alpha_i;
        A(ii,ii)   = -(alpha_i+gamma_i);
        A(ii,ii+1) = gamma_i;
    end
end

% Create vector b
b=zeros(N+2,1);
for ii=1:N
    if(ii==1)
        b(ii)   = u0;
    elseif(ii==N+2)
        b(ii)   = (f(x(ii)) + f(x(ii)+dx))/2.0 + 2*uN;
    else
        %        b(ii)  =   f((x(ii) + x(ii+1))/2.0);       % Midpoint rule
        b(ii)   =   (f(x(ii)) + f(x(ii+1)))/2.0;            % Trepozoidal rule
    end
end

% u=zeros(N,1);
u = A\b;