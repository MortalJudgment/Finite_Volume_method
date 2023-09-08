% Error in H1
function errorH1 = normH1(w_dis,w_ex,x,y,x_cp,y_cp,Width_x,Width_y)
Nx = Width_x;
Ny = Width_y;
e_ij = zeros(Nx+2,Ny+2);
for jj=1:Ny+2
    for ii=1:Nx+2
        e_ij(ii,jj) = abs(w_dis(ii,jj) - w_ex(ii,jj));
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
        e2 = e2 + (x(ii)-x(ii-1))/(y_cp(jj+1)-y_cp(jj))*(e_ij(ii,jj+1)-e_ij(ii,jj))^2;
    end
end
errorH1 = e1+e2;
errorH1 = sqrt(errorH1);