function r = Aj(x,x_cp,y,y_cp,Nx,jj)
r = zeros(Nx+2,Nx+2);
for ii=1:Nx+2
    if ii==1
        r(ii,ii) = 1;
    elseif ii==Nx+2
        r(ii,ii) = 1;
    else
        a_i     = 1.0/((x(ii)-x(ii-1))*(x_cp(ii)-x_cp(ii-1)));
        b_i     = 1.0/((x(ii)-x(ii-1))*(x_cp(ii+1)-x_cp(ii)));
        c_j     = 1.0/((y(jj)-y(jj-1))*(y_cp(jj)-y_cp(jj-1)));
        d_j     = 1.0/((y(jj)-y(jj-1))*(y_cp(jj+1)-y_cp(jj)));
        s_ij    = a_i + b_i + c_j + d_j;
        
        r(ii,ii-1) = -a_i;
        r(ii,ii) = s_ij;
        r(ii,ii+1) = -b_i;
    end
end