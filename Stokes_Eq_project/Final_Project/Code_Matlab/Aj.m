function r = Aj(x,x_cp,y,y_cp,Nx,jj)
r = zeros(Nx,Nx);
for ii=1:Nx
    a_i     = 1.0/((x(ii+1)-x(ii))*(x_cp(ii+1)-x_cp(ii)));
    b_i     = 1.0/((x(ii+1)-x(ii))*(x_cp(ii+2)-x_cp(ii+1)));
    c_j     = 1.0/((y(jj+1)-y(jj))*(y_cp(jj+1)-y_cp(jj)));
    d_j     = 1.0/((y(jj+1)-y(jj))*(y_cp(jj+2)-y_cp(jj+1)));
    s_ij    = a_i + b_i + c_j + d_j;
    
    if ii==1
        r(ii,ii) = s_ij;
        r(ii,ii+1) = -b_i;
    elseif ii==Nx
        r(ii,ii-1) = -a_i;
        r(ii,ii) = s_ij;
    else
        r(ii,ii-1) = -a_i;
        r(ii,ii) = s_ij;
        r(ii,ii+1) = -b_i;
    end
end