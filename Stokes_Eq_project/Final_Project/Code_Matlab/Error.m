%%
% Order of convergence
% Finite Volume Method for Stokes problem
%%
clc
clear all
close all

type = 1;
%----------------%
% Domain of x-axis
ax = 0.0;
bx = 1.0;
% Domain of y-axis
ay = 0.0;
by = 1.0;
%----------------%
%%
Nx = 5;
Ny = Nx;
mtis = 4;
ll = zeros(mtis,1);
% error of u
errorL2_u = zeros(mtis,1); % error in L2
errorH1_u = zeros(mtis,1); % error in H1
% error of v
errorL2_v = zeros(mtis,1); % error in L2
errorH1_v = zeros(mtis,1); % error in H1
% error of p
errorL2_p = zeros(mtis,1); % error in L2
errorH1_p = zeros(mtis,1); % error in H1
%%
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
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------- Creare the Matrix --------------------------%
    % Matrix present Divergence operator
    A=zeros(Nx*Ny,Nx*Ny);
    for jj = 1:Ny
        if jj==1
            A((jj-1)*Nx+1:jj*Nx,(jj-1)*Nx+1:(jj+1)*Nx) = [Aj(x,x_cp,y,y_cp,Nx,jj),Dj(y,y_cp,Nx,jj)];
        elseif jj==Ny
            A((jj-1)*Nx+1:jj*Nx,(jj-2)*Nx+1:jj*Nx) = [Cj(y,y_cp,Nx,jj),Aj(x,x_cp,y,y_cp,Nx,jj)];
        else
            A((jj-1)*Nx+1:jj*Nx,(jj-2)*Nx+1:(jj+1)*Nx) = [Cj(y,y_cp,Nx,jj),Aj(x,x_cp,y,y_cp,Nx,jj),Dj(y,y_cp,Nx,jj)];
        end
    end
    %%
    % Gradient matrix operator on x
    B1=zeros(Nx*Ny,Nx*Ny);
    for jj = 1:Ny
        B1((jj-1)*Nx+1:jj*Nx,(jj-1)*Nx+1:jj*Nx) = E1j(x,Nx,jj);
    end
    %%
    % Gradient matrix operator on y
    B2=zeros(Nx*Ny,Nx*Ny);
    for jj = 1:Ny
        if jj==1
            B2((jj-1)*Nx+1:jj*Nx,jj*Nx+1:(jj+1)*Nx) = E2j(y,Nx,jj);
            B2((jj-1)*Nx+1:jj*Nx,(jj-1)*Nx+1:jj*Nx) = -E2j(y,Nx,jj);
        elseif jj==Ny
            B2((jj-1)*Nx+1:jj*Nx,(jj-1)*Nx+1:jj*Nx) = E2j(y,Nx,jj);
            B2((jj-1)*Nx+1:jj*Nx,(jj-2)*Nx+1:(jj-1)*Nx) = -E2j(y,Nx,jj);
        else
            B2((jj-1)*Nx+1:jj*Nx,(jj-2)*Nx+1:(jj-1)*Nx) = -E2j(y,Nx,jj);
            B2((jj-1)*Nx+1:jj*Nx,jj*Nx+1:(jj+1)*Nx) = E2j(y,Nx,jj);
        end
    end
    %%
    %--------------------------------%
    % Value of function f(x,y)
    % Value of function f(x,y)
    F1 = zeros(Nx*Ny,1);
    F2 = zeros(Nx*Ny,1);
    for jj=1:Ny
        for ii=1:Nx
            F1((jj-1)*Ny+ii) = (f(0.5*x(ii)+0.5*x(ii+1),y(jj),type)+f(x(ii+1),0.5*y(jj)+0.5*y(jj+1),type)+f(x(ii),0.5*y(jj)+0.5*y(jj+1),type)+f(0.5*x(ii)+0.5*x(ii+1),y(jj+1),type))/4.0;
            F2((jj-1)*Ny+ii) = (g(0.5*x(ii)+0.5*x(ii+1),y(jj),type)+g(x(ii+1),0.5*y(jj)+0.5*y(jj+1),type)+g(x(ii),0.5*y(jj)+0.5*y(jj+1),type)+g(0.5*x(ii)+0.5*x(ii+1),y(jj+1),type))/4;
        end
    end
    %%
    % Boundary conditions
    % u, v = 0;
    
    %% Type of system equation
    %
    % / A  0  B1 \ / u \  _  /  F1  \
    % \ 0  A  B2 / \ v /  _  \  F2  /
    %  \B1 B2 0 /   \p/       \  0 /
    % by let AA = [A,0;0,A]
    %         B = [B1;B2]
    %        x1 = [u;v]
    %        x2 = [p]
    %         F = [F1;F2]
    % we have new equation had typed:
    % /AA  B \ / x1 \     /F\
    % \ B* 0 / \ x2 /  =  \0/
    AA = [A,0.*A;0.*A,A];
    B = [B1;B2];
    FF = [F1;F2];
    
    %-------  Using Uzawa iteration solving the system equation upper -------%
    [x1,x2] = Uzawa_iteration(AA,B,FF,Nx,Ny,dx);
    u = x1(1:Nx*Ny);
    v = x1(Nx*Ny+1:2*Nx*Ny);
    p = x2;
    
    %%
    %---------- Dicrete solution -----------%
    %-- u dicrete --%
    u1 = zeros(Nx,Ny);
    kk = 1;
    for jj=1:Ny
        for ii=1:Nx
            u1(ii,jj) = u(kk);
            kk = kk + 1;
        end
    end
    u_dis = zeros(Nx+2,Ny+2);
    u_dis(2:Nx+1,2:Ny+1) = u1;
    
    %-- v dicrete --%
    v1 = zeros(Nx,Ny);
    kk = 1;
    for jj=1:Ny
        for ii=1:Nx
            v1(ii,jj) = v(kk);
            kk = kk + 1;
        end
    end
    v_dis = zeros(Nx+2,Ny+2);
    v_dis(2:Nx+1,2:Ny+1) = v1;
    
    %--- p dicrete ---%
    p1 = zeros(Nx,Ny);
    kk = 1;
    for jj=1:Ny
        for ii=1:Nx
            p1(ii,jj) = p(kk);
            kk = kk + 1;
        end
    end
    p_dis = p1;
    %%
    %------------------------------------------------%
    %---------- Exact solution -----------%
    %--- u exact ---%
    u_exact = zeros(Nx+2,Ny+2,2);
    for jj=1:Ny+2
        for ii=1:Nx+2
            u_exact(ii,jj,:) = uex(x_cp(ii),y_cp(jj),type);
        end
    end
    u_ex = u_exact(:,:,1);
    v_ex = u_exact(:,:,2);
    %--- p exact ---%
    p_ex = zeros(Nx,Ny);
    for jj=1:Ny
        for ii=1:Nx
            p_ex(ii,jj) = pex(x_cp(ii+1),y_cp(jj+1),type);
        end
    end
    %%
    %---------------- Error in L2 ----------------%
    
    errorL2_u(iter) = normL2(u_dis,u_ex,x,y,x_cp,y_cp,Nx,Ny);
    errorL2_v(iter) = normL2(v_dis,v_ex,x,y,x_cp,y_cp,Nx,Ny);
    errorL2_p(iter) = normL2(p_dis,p_ex,x(2:Nx-1),y(2:Ny-1),x_cp(2:Nx+1),y_cp(2:Nx+1),Nx-3,Ny-3);
    %%
    %---------------- Error in H1 ----------------%
    
    errorH1_u(iter) = normH1(u_dis,u_ex,x,y,x_cp,y_cp,Nx,Ny);
    errorH1_v(iter) = normH1(v_dis,v_ex,x,y,x_cp,y_cp,Nx,Ny);
    errorH1_p(iter) = normH1(p_dis,p_ex,x(2:Nx-1),y(2:Ny-1),x_cp(2:Nx+1),y_cp(2:Nx+1),Nx-3,Ny-3);
    %%
    ll(iter) = (Nx+Ny)/2.0;
    Nx = Nx*2;
    Ny = Nx;
end
figure
plot(log(ll),-log(errorL2_u),'r',log(ll)-1,-log(errorH1_u)-2,'y',log(ll),1*log(ll)-1.5,'b',log(ll),2*log(ll)-2,'g')
legend('normL2','normH1','1x','2x')
title('Order of approximate for u')
figure
plot(log(ll),-log(errorL2_v),'r',log(ll)-1,-log(errorH1_v)-2,'y',log(ll),1*log(ll)-1.5,'b',log(ll),2*log(ll)-2,'g')
legend('normL2','normH1','1x','2x')
title('Order of approximate for v')
figure
plot(log(ll),-log(errorL2_p)+2,'r',log(ll),-log(errorH1_p),'y',log(ll),1*log(ll)-1.5,'b',log(ll),2*log(ll)-2,'g')
legend('normL2','normH1','1x','2x')
title('Order of approximate for p')