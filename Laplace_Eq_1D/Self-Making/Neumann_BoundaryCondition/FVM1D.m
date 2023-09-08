%-------------------------------------------------------------------------%
%----- Solve 1D Laplace equation: -----% 
%                 -u_xx = f(x)          ,in [a,b]
%-------------------------------------------------------------------------%
clc
clear all
close all
%----------------%
% Domain of x
ax = 0.0;
bx = 1.0;
% Boundary Condition
u0 = 0;
uN = 0;
%----------------%
N = 4;                      % Number of control volume
M = 6;                      % number of iteration when refine mesh
%----------------%
norml2 = zeros(M,1);        % norm in L_2
normh1 = zeros(M,1);        % norrm in H_1

ll = zeros(M,1);            % vector contain control value


for jj = 1:M
    %---------- Dicrete solution -----------%
    u_dis = FVM1D(ax,bx,N,u0,uN);
    
    %----------- Exact solution ------------%
    u_ex = zeros(N+2,1);
    for ii=1:N+2
        u_ex(ii) = u_exact(x_cp(ii));
    end   
    %---------------------------------------%
    %-------------- Drawing ----------------%
    %---------------------------------------%
    figure
    plot(x_cp,u_dis,'red',x_cp,u_ex);
    
    % Error norm in L2
    for ii=1:N
        norml2(jj) = norml2(jj)+(u_dis(ii+1)-u_ex(ii+1))^2*(x(ii+1)-x(ii));
    end
    norml2(jj)=sqrt(norml2(jj));
    % Error norm in H1
    for ii=1:N+1
        normh1(jj) = normh1(jj)+((u_dis(ii+1)-u_ex(ii+1))-(u_dis(ii)-u_ex(ii)))^2/(x_cp(ii+1)-x_cp(ii));
    end
    normh1(jj) = sqrt(normh1(jj));
    
    % Update ll(jj)
    ll(jj) = N;
    % Update number of control volume
    N=2*N;
end

figure
plot(log(ll),-log(norml2),'r', log(ll), -log(normh1),'blue', log(ll),1.5*log(ll)+2, 'black', log(ll), 2*log(ll)+1.5,'green',  log(ll), log(ll),'y');
title('Error');
legend('L^2 Norm', 'H^1 norm', '3/2x', '2x')