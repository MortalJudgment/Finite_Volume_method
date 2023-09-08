function r = Dj(y,y_cp,Nx,jj)
r = zeros(Nx,Nx);
for ii=1:Nx+2
    d_j     = 1.0/((y(jj)-y(jj-1))*(y_cp(jj+1)-y_cp(jj)));
    r(ii,ii) = -d_j;
end