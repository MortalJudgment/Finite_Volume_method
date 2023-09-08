function r = Dj(y,y_cp,Nx,jj)
r = zeros(Nx,Nx);
for ii=1:Nx
    d_j     = 1.0/((y(jj+1)-y(jj))*(y_cp(jj+2)-y_cp(jj+1)));
    r(ii,ii) = -d_j;
end