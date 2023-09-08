%% Solve equation -lapla(u)=f(x) with the Dirichlet boundary condition 
clear all
close all
clc
%% Initial informations
ax=0.0;
bx=1.0;
cases=3;
N=5; %number of mesh points of first mesh
number_mesh=4;
number_mesh_point=zeros(number_mesh,1);
norm_max=zeros(number_mesh,1);
norm_l2=zeros(number_mesh,1);
norm_maxh1=zeros(number_mesh,1);
norm_h1=zeros(number_mesh,1);
%% Solve discrite solution and refine mesh

for inumber_mesh=1:number_mesh
    
    number_mesh_point(inumber_mesh)=N;
    h=1/N;
%% Create mesh point    
    x=linspace(ax,bx,N+1)';
    y=x;
% % Create Control Point
    z=zeros(N+2,1);
    z(N+2)=bx;
    for i=2:N+1
         z(i)=0.5*x(i-1)+0.5*x(i);
    end
    t=z;

 %% Create matrix A   
    A=zeros((N)^2,(N)^2);
    a=zeros(N,1);
    b=zeros(N,1);
    c=zeros(N,1);
    d=zeros(N,1);
    s=zeros(N);
    for i=1:N
        a(i)=1/((x(i+1)-x(i))*(z(i+1)-z(i)));
        b(i)=1/((x(i+1)-x(i))*(z(i+2)-z(i+1)));
        c(i)=1/((y(i+1)-y(i))*(t(i+1)-t(i)));
        d(i)=1/((y(i+1)-y(i))*(t(i+2)-t(i+1)));
    end
    for i=1:N
        for j=1:N
            s(i,j)=a(i)+b(i)+c(j)+d(j);
        end
    end
    B=zeros(N,N,N);
    C=zeros(N,N,N-1);
    D=zeros(N,N,N-1);
    for k=1:N
        for i=1:N
            if(i==1)
                    B(i,i,k)=s(i,k);
                    B(i,i+1,k)=-b(i);
                elseif(i==N)
                    B(i,i,k)=s(i,k);
                    B(i,i-1,k)=-a(i);
                else B(i,i+1,k)=-b(i);
                     B(i,i-1,k)=-a(i);
                     B(i,i,k)=s(i,k);
            end
        end
    end

    for k=1:N-1
        for i=1:N
            C(i,i,k)=-c(k+1);
        end
    end
    
    for k=1:N-1
        for i=1:N
            D(i,i,k)=-d(k);
        end
    end
    
    for i=1:N:N^2
        if(i==1)
        A(i:(i+N-1),i:(i+N-1))=B(:,:,i);
        A(i:(i+N-1),(i+N):(i+N+N-1))=D(:,:,i);
        elseif(i-1+N==N^2)
        A(i:(i+N-1),i:(i+N-1))=B(:,:,(i-1)/N+1);
        A(i:(i+N-1),(i-N):i-1)=C(:,:,(i-1)/N+1-1);
        else
        A(i:(i+N-1),i:(i+N-1))=B(:,:,(i-1)/N+1);
        A(i:(i+N-1),(i+N):(i+N+N-1))=D(:,:,(i-1)/N+1);
        A(i:(i+N-1),(i-N):i-1)=C(:,:,(i-1)/N+1-1);
        end
    end
    
    K=zeros(3*(N)^2,3*(N)^2);
    K(1:N^2,1:N^2)=A;
    K(N^2+1:2*N^2,N^2+1:2*N^2)=A;
    
    A=zeros((N)^2,(N)^2);
    B=zeros(N,N,N);
    for k=1:N
        for i=1:N
            if(i==1)
                    B(i,i,k)=-1;
                    B(i,i+1,k)=1;
                elseif(i==N)
                    B(i,i,k)=1;
                    B(i,i-1,k)=-1;
            else
                B(i,i+1,k)=1;
                B(i,i-1,k)=-1;
                B(i,i,k)=0;
            end
        end
    end
    
    B=(1/(2*h))*B;
    for i=1:N:N^2
        if(i==1)
        A(i:(i+N-1),i:(i+N-1))=B(:,:,i);
        elseif(i-1+N==N^2)
        A(i:(i+N-1),i:(i+N-1))=B(:,:,(i-1)/N+1);
        else
        A(i:(i+N-1),i:(i+N-1))=B(:,:,(i-1)/N+1);
        end
    end
    
    K(1:N^2,2*N^2+1:3*N^2)=A;
     
    
    A=zeros((N)^2,(N)^2);
    
    for i=1:N:N^2
        if(i==1)
        A(i:(i+N-1),i:(i+N-1))=-eye(N);
        A(i:(i+N-1),(i+N):(i+N+N-1))=eye(N);
        elseif(i-1+N==N^2)
        A(i:(i+N-1),i:(i+N-1))=eye(N);
        A(i:(i+N-1),(i-N):i-1)=-eye(N);
        else
        A(i:(i+N-1),i:(i+N-1))=0;
        A(i:(i+N-1),(i+N):(i+N+N-1))=eye(N);
        A(i:(i+N-1),(i-N):i-1)=-eye(N);
        end
    end
    
    K(N^2+1:2*N^2,2*N^2+1:3*N^2)=(1/(2*h))*A;
    
    A=zeros((N)^2,(N)^2);
    B=zeros(N,N,N);
    for k=1:N
        for i=1:N
            if(i==1)
                   B(i,i,k)=1;
                   B(i,i+1,k)=1;
            elseif(i==N)
                   B(i,i,k)=-1;
                   B(i,i-1,k)=-1;
            else
                B(i,i+1,k)=1;
                B(i,i-1,k)=-1;
                B(i,i,k)=0;
            end
        end
    end
    
    for i=1:N:N^2
        if(i==1)
        A(i:(i+N-1),i:(i+N-1))=B(:,:,i);
        elseif(i-1+N==N^2)
        A(i:(i+N-1),i:(i+N-1))=B(:,:,(i-1)/N+1);
        else
        A(i:(i+N-1),i:(i+N-1))=B(:,:,(i-1)/N+1);
        end
    end
    
    K(2*N^2+1:3*N^2,1:N^2)=-(1/(2*h))*A;

    
    A=zeros((N)^2,(N)^2);
    
    for i=1:N:N^2
        if(i==1)
        A(i:(i+N-1),i:(i+N-1))=eye(N);
        A(i:(i+N-1),(i+N):(i+N+N-1))=eye(N);
        elseif(i-1+N==N^2)
        A(i:(i+N-1),i:(i+N-1))=-eye(N);
        A(i:(i+N-1),(i-N):i-1)=-eye(N);
        else
        A(i:(i+N-1),(i-N):i-1)=-eye(N);
        A(i:(i+N-1),(i+N):(i+N+N-1))=eye(N);
        end
    end
    
    K(2*N^2+1:3*N^2,N^2+1:2*N^2)=-(1/(2*h))*A;
    
%% Create vector b    
    f=zeros((N)*(N),1);
    g=zeros((N)*(N),1);
    ff=zeros(3*(N)*(N),1);
    for j=1:N
        for i=1:N
              f((j-1)*N+i)=(function_f1(0.5*x(i)+0.5*x(i+1),y(j),cases)+function_f1(x(i+1),0.5*y(j)+0.5*y(j+1),cases)+function_f1(x(i),0.5*y(j)+0.5*y(j+1),cases)+function_f1(0.5*x(i)+0.5*x(i+1),y(j+1),cases))/4;
        end
    end
    
    for j=1:N
        for i=1:N
              g((j-1)*N+i)=(function_f2(0.5*x(i)+0.5*x(i+1),y(j),cases)+function_f2(x(i+1),0.5*y(j)+0.5*y(j+1),cases)+function_f2(x(i),0.5*y(j)+0.5*y(j+1),cases)+function_f2(0.5*x(i)+0.5*x(i+1),y(j+1),cases))/4;
        end
    end
    
    ff(1:N^2)=f;
    ff(N^2+1:2*N^2)=g;
%% Solve discrete solution


%     K(3*N^2,2*N^2+1:3*N^2)=-1;
%     u=K\ff;



%% Solve Ax=b by Uzawa interation by compute completely x_2 then compute x_1 

    R=K(1:2*N^2,2*N^2+1:3*N^2);
    T=K(1:2*N^2,1:2*N^2);
    S=(R')*(T^-1)*(R);
    SS=(T^-1)*(R);
    b_1=ff(1:2*N^2,1);
    b_vector=R'*(T^-1)*b_1;
    u_guess(1:N^2,1)=1;
    x_1=T^-1*(b_1-R*u_guess);
    p=conjgrad(S,b_vector,u_guess);
     
    sum=0.0;
    for i=1:N^2
         sum=sum+p(i)*h^2;
    end
    for i=1:N^2
        p(i)=p(i)-sum;
    end
        
    u1=T^-1*(b_1-R*p);
    u=[u1;p];
    
%% Solve Ax=b by Uzawa interation by compute  x_2 alongside x_1 

    
    %% Create discrete solution with boundary 
    u_dis=zeros((N+2)^2,1);
    ii=1;
    for j=1:N+2
        for i=1:N+2
           if(i==1 || i==N+2 || j==1 || j==N+2)
               u_dis((j-1)*(N+2)+i)=0;
           else
               u_dis((j-1)*(N+2)+i)=u(ii);
               ii=ii+1;
           end
        end
    end
    
    v_dis=zeros((N+2)^2,1);
    for j=1:N+2
        for i=1:N+2
           if(i==1 || i==N+2 || j==1 || j==N+2)
               v_dis((j-1)*(N+2)+i)=0;
           else
               v_dis((j-1)*(N+2)+i)=u(ii);
               ii=ii+1;
           end
        end
    end
    
    p_dis=zeros((N)^2,1);
    for j=1:N
        for i=1:N
            p_dis((j-1)*(N)+i)=u(ii);
            ii=ii+1;
        end
    end
        
%% Get exact solution    
    u_ex=zeros((N+2)^2,1);
    for j=1:N+2
        for i=1:N+2
            u_ex((j-1)*(N+2)+i)=exact_solution_u(z(i),t(j),cases);
        end
    end
    
    v_ex=zeros((N+2)^2,1);
    for j=1:N+2
        for i=1:N+2
            v_ex((j-1)*(N+2)+i)=exact_solution_v(z(i),t(j),cases);
        end
    end
    
    p_ex=zeros(N^2,1);
    for j=2:N+1
        for i=2:N+1
            p_ex((j-2)*(N)+i-1)=exact_pressure(z(i),t(j),cases);
        end
    end
    
%% Figure exact and dicrete solutions    
%     figure
    [X,Y]=meshgrid(z,t);
    U_dis=zeros(size(X));
    U_ex=zeros(size(X));
    
    for i=1:N+2
        U_dis(i,1:N+2)=u_dis((i-1)*(N+2)+1:(i-1)*(N+2)+1+N+1,1);
        U_ex(i,1:N+2)=u_ex((i-1)*(N+2)+1:(i-1)*(N+2)+1+N+1,1);
    end
    
    V_dis=zeros(size(X));
    V_ex=zeros(size(X));
    for i=1:N+2
        V_dis(i,1:N+2)=v_dis((i-1)*(N+2)+1:(i-1)*(N+2)+1+N+1,1);
        V_ex(i,1:N+2)=v_ex((i-1)*(N+2)+1:(i-1)*(N+2)+1+N+1,1);
    end
    
    P_dis=zeros(N);
    P_ex=zeros(N);
    [XX,YY]=meshgrid(z(2:N+1),t(2:N+1));
    for i=1:N
        P_dis(i,1:N)=p_dis((i-1)*(N)+1:(i-1)*(N)+N,1);
        P_ex(i,1:N)=p_ex((i-1)*(N)+1:(i-1)*(N)+N,1);
    end
    
%     figure
%     subplot(1,2,1);
%     surf(Y,X,U_ex);
% %     view(0,90)
% %     colorbar
%     xlabel('x');ylabel('y');zlabel('value');
%     title('Exact solution u');
%     
%     subplot(1,2,2);
%     surf(Y,X,U_dis);
% %     view(0,90)
% %     colorbar
%     title('Discrete solution u');
%     xlabel('x');ylabel('y');zlabel('value');
%     
%     figure
%     subplot(1,2,1);
%     surf(Y,X,V_ex);
% %     view(0,90)
% %     colorbar
%     xlabel('x');ylabel('y');zlabel('value');
%     title('Exact solution v');
%     subplot(1,2,2);
%     surf(Y,X,V_dis);
% %     view(0,90)
% %     colorbar
%     title('Discrete solution v');
%     xlabel('x');ylabel('y');zlabel('value');
%     
%     figure
%     subplot(1,2,1);
%     surf(YY,XX,P_ex);
% %     contourf(YY,XX,P_ex);
% %     view(0,90)
% %     colorbar
%     xlabel('x');ylabel('y');zlabel('value');
%     title('Exact solution p');
%     subplot(1,2,2);
%     surf(YY,XX,P_dis);
% %     contourf(YY,XX,P_dis);
% %     view(0,90)
% %     colorbar
%     title('Discrete solution p');
%     xlabel('x');ylabel('y');zlabel('value');
%     
% 
%     
  
%%  Calculate the error on L^2 

    norm_l2_u(inumber_mesh)=0;
    for i=2:N+1
        for j=2:N+1
           norm_l2_u(inumber_mesh)=norm_l2_u(inumber_mesh)+(U_dis(i,j)-U_ex(i,j))^2*(x(i)-x(i-1))*(y(j)-y(j-1)); 
        end
    end
    norm_l2_u(inumber_mesh)=(norm_l2_u(inumber_mesh))^(1/2)

    norm_l2_v(inumber_mesh)=0;
    for i=2:N+1
        for j=2:N+1
           norm_l2_v(inumber_mesh)=norm_l2_v(inumber_mesh)+(V_dis(i,j)-V_ex(i,j))^2*(x(i)-x(i-1))*(y(j)-y(j-1)); 
        end
    end
    norm_l2_v(inumber_mesh)=(norm_l2_v(inumber_mesh))^(1/2)

    norm_l2_p(inumber_mesh)=0;
    
    for i=1:N
        for j=1:N
           norm_l2_p(inumber_mesh)=norm_l2_p(inumber_mesh)+(P_dis(i,j)-P_ex(i,j))^2*(h^2); 
        end
    end
    norm_l2_p(inumber_mesh)=(norm_l2_p(inumber_mesh))^(1/2)
    
    %% Calculate the error on H1

    norm_h1_u(inumber_mesh)=0;
    for i=1:N
        for j=2:N+2
        norm_h1_u(inumber_mesh)=norm_h1_u(inumber_mesh)+(((((U_dis(i+1,j)-U_ex(i+1,j))-(U_dis(i,j)-U_ex(i,j))))^2)/(z(j)-z(j-1)))*(y(i+1)-y(i));
        end
    end
    for j=1:N+1
        for i=2:N+1
        norm_h1_u(inumber_mesh)=norm_h1_u(inumber_mesh)+(((((U_dis(i,j+1)-U_ex(i,j+1))-(U_dis(i,j)-U_ex(i,j))))^2)/(t(j+1)-t(j)))*(x(i)-x(i-1));
        end
    end
    norm_h1_u(inumber_mesh)=(norm_h1_u(inumber_mesh))^(1/2)

    norm_h1_v(inumber_mesh)=0;
    for i=1:N
        for j=2:N+2
        norm_h1_v(inumber_mesh)=norm_h1_v(inumber_mesh)+(((((V_dis(i+1,j)-V_ex(i+1,j))-(V_dis(i,j)-V_ex(i,j))))^2)/(z(j)-z(j-1)))*(y(i+1)-y(i));
        end
    end
    
    for j=1:N+1
        for i=2:N+1
        norm_h1_v(inumber_mesh)=norm_h1_v(inumber_mesh)+(((((V_dis(i,j+1)-V_ex(i,j+1))-(V_dis(i,j)-V_ex(i,j))))^2)/(t(j+1)-t(j)))*(x(i)-x(i-1));
        end
    end
    norm_h1_v(inumber_mesh)=(norm_h1_v(inumber_mesh))^(1/2)
     
%% Refine mesh (increse mesh point)    
    N=2*N;
end

%% Figure for errors respect to number of mesh point
figure
plot( log(number_mesh_point), -log(norm_l2_u), 'yellow',...
     log(number_mesh_point), -log(norm_h1_u)-1, 'magenta',...
    log(number_mesh_point), 2*log(number_mesh_point)+1,'red',log(number_mesh_point), 1.5*log(number_mesh_point)-1,'blue');
% xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Errors on u');
legend('normL2 u','normH1 u','2x','1.5x'); 
figure
plot( log(number_mesh_point), -log(norm_l2_v), 'yellow',...
     log(number_mesh_point), -log(norm_h1_v)-0.5, 'magenta',...
    log(number_mesh_point), 2*log(number_mesh_point)+1,'red',log(number_mesh_point), 1.5*log(number_mesh_point)-1.5,'blue');
% xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Errors on v');
legend('normL2 v','normH1 v','2x','1.5x'); 
figure

plot( log(number_mesh_point), -log(norm_l2_p), 'yellow',...
    log(number_mesh_point), 2*log(number_mesh_point),'black',log(number_mesh_point), 1.5*log(number_mesh_point)-1.5,'blue',...
    log(number_mesh_point), log(number_mesh_point)+1,'green');
% xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Errors on p');
legend('normL2 p','2x','1.5x','x'); 


