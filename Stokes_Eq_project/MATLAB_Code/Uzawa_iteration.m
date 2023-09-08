%----------  Using Uzawa iteration solving the system equation  ----------%
% / A  B \ / x1 \     /F\
% \ B* 0 / \ x2 /  =  \0/

function [x1,x2] = Uzawa_iteration(A,B,F,Nx,Ny,dx)
% Schur complement
S = (B')*A^(-1)*B;
% S is symmetric positive-definete we can apply standard iterative methods
% like the conjugate gradient method to solve equation
% S.x2 = B*(A^(-1))*F to compute x2

% Conjugate gradient method:
% Initial guess
x2 = eye(Nx*Ny,1);
q = B'*(A^-1)*F;
% Start a lope
r1 = q - S * x2;
p2 = r1;
rs = r1' * r1;

for i = 1:length(q)
    a2 = S * p2;
    alpha = rs / (p2' * a2);
    x2 = x2 + alpha * p2;
    r1 = r1 - alpha * a2;
    r2 = r1' * r1;
    if sqrt(r2) < 1e-10
        break;
    end
    p2 =  r1 + (r2 / rs) * p2;
    rs = r2;
end

sum=0.0;
for i=1:Nx*Ny
    sum=sum+x2(i)*dx^2;
end
for i=1:Nx*Ny
    x2(i)=x2(i)-sum;
end
x1 = A^-1*(F-B*x2);
end