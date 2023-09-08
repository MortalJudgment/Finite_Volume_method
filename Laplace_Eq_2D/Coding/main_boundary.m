%-------------------------------------------------------------------------%
%----- Solve 1D Laplace equation: -----%
%                 -u_xx = f(x)          ,in [a,b]
%-------------------------------------------------------------------------%
clc
clear all
close all
%----------------%
% Domain of x-axis
ax = 0.0;
bx = 1.0;
% Domain of y-axis
ay = 0.0;
by = 1.0;
%----------------%
Nx = 5;
Ny = Nx;
mtis = 4;
ll = zeros(mtis,1);
errorL2 = zeros(mtis,1); % error in L2
errorH1 = zeros(mtis,1); % error in H1
for iter=1:mtis
    dx = (bx-ax)/Nx;
    dy = (by-ay)/Ny;
    % Create the mesh point
    x = zeros(Nx+1,1);
    y = zeros(Ny+1,1);
    for ii=1:Nx+1
        x(ii) = ax+(ii-1)*dx;
    end
    for ii=1:Ny+1
        y(ii) = ay+(ii-1)*dy;
    end
    % Create control point
    x_cp = zeros(Nx+2,1);
    y_cp = zeros(Ny+2,1);
    for ii=1:Nx+2
        if(ii==1)
            x_cp(ii) = x(ii);
        elseif(ii==Nx+2)
            x_cp(ii) = x(ii-1);
        else
            x_cp(ii) = (1*x(ii-1)+1*x(ii))/2.0;
        end
    end
    for ii=1:Ny+2
        if(ii==1)
            y_cp(ii) = y(ii);
        elseif(ii==Ny+2)
            y_cp(ii) = y(ii-1);
        else
            y_cp(ii) = (1*y(ii-1)+1*y(ii))/2.0;
        end
    end
    % Creare the Matrix
    A=zeros(Nx*Ny,Nx*Ny);
    for jj = 1:Ny
        if jj==1
            A((jj-1)*Nx+1:jj*Nx,(jj-1)*Nx+1:(jj+1)*Nx) = [Aj(x,x_cp,y,y_cp,Nx,jj),Dj(y,y_cp,Nx,jj)];
        elseif jj==Ny
            A((jj-1)*Nx+1:jj*Nx,(jj-2)*Nx+1:(jj)*Nx) = [Cj(y,y_cp,Nx,jj),Aj(x,x_cp,y,y_cp,Nx,jj)];
        else
            A((jj-1)*Nx+1:jj*Nx,(jj-2)*Nx+1:(jj+1)*Nx) = [Cj(y,y_cp,Nx,jj),Aj(x,x_cp,y,y_cp,Nx,jj),Dj(y,y_cp,Nx,jj)];
        end
    end
    
    b=zeros(Nx*Ny,1);
    kk=1;
    for jj=1:Ny
        for ii=1:Nx
            b(kk)   =   (f(x(ii),y(jj)) + f(x(ii+1),y(jj)) + f(x(ii),y(jj+1)) + f(x(ii+1),y(jj+1)))/4.0;
            kk = kk + 1;
        end
    end
    u=A\b;
%----------- Exact solution ------------%
    u_ex = zeros((Nx+2),(Ny+2));
    for jj=1:Ny+2
        for ii=1:Nx+2
            u_ex(ii,jj) = u_exact(x_cp(ii),y_cp(jj));
        end
    end  
%---------- Dicrete solution -----------%
    u1 = zeros(Nx,Ny);
    kk = 1;
    for jj=1:Ny
        for ii=1:Nx
            u1(ii,jj) = u(kk);
            kk = kk + 1;
        end
    end
    u_dis = zeros((Nx+2),(Ny+2));
    u_dis(2:Nx+1,2:Ny+1) = u1;
%---------------------------------------%
%-------------- Drawing ----------------%
%---------------------------------------%
    figure
    subplot(1,2,1)
    h1 = surf(x_cp,y_cp,u_dis);
    title('Discrete solution')
    subplot(1,2,2)
    h2 = surf(x_cp,y_cp,u_ex);
    title('Exactly solution')
%---------------- Error in L2 ----------------%
    e_ij = zeros(Nx+2,Ny+2);
    for jj=1:Ny+2
        for ii=1:Nx+2
            e_ij(ii,jj) = abs(u_dis(ii,jj) - u_ex(ii,jj));
        end
    end
    e1 = 0;
    for jj=2:Ny+1
        for ii=1:Nx+1
            e1 = e1 + (y(jj)-y(jj-1))/(x_cp(ii+1)-x_cp(ii))*(e_ij(ii+1,jj)-e_ij(ii,jj))^2;
        end
    end
    e2 = 0;
    for jj=1:Ny+1
        for ii=2:Nx+1
            e2 = e2 + (x(ii)-x(ii-1))/(y_cp(ii+1)-y_cp(ii))*(e_ij(ii,jj+1)-e_ij(ii,jj))^2;
        end
    end
    errorL2(iter) = e1+e2;
    errorL2(iter) = sqrt(errorL2(iter));
%---------------- Error in H1 ----------------%
    errorH1(iter) = 0;
    for jj=2:Ny+1
        for ii=2:Nx+1
            errorH1 = errorH1 + (e_ij(ii,jj))^2*(x(ii)-x(ii-1))*(y(ii)-y(ii-1));
        end
    end
    errorH1(iter) = sqrt(errorH1(iter));
    ll(iter) = (Nx+Ny)/2.0;
    Nx = Nx*2;
    Ny = Nx;
end
figure
plot(log(ll),-log(errorL2),'r',log(ll),1.5*log(ll)+4,'b',log(ll),-log(errorH1),'y',log(ll),2*log(ll)+4,'g')
legend('normL2','1.5x','normH1','2x')
title('Error')