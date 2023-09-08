% Error in L2
function errorL2 = normL2(w_dis,w_ex,x,y,x_cp,y_cp,Width_x,Width_y)
Nx = Width_x;
Ny = Width_y;
errorL2 = 0;
e_ij = zeros(Nx+2,Ny+2);
for jj=1:Ny+2
    for ii=1:Nx+2
        e_ij(ii,jj) = (w_dis(ii,jj) - w_ex(ii,jj));
    end
end
for jj=2:Ny+1
    for ii=2:Nx+1
        errorL2 = errorL2 + (e_ij(ii,jj))^2*(x(ii)-x(ii-1))*(y(ii)-y(ii-1));
    end
end
errorL2 = sqrt(errorL2);