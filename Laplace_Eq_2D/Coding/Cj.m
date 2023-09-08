function r = Cj(y,y_cp,Nx,jj)
r = zeros(Nx,Nx);
for ii=1:Nx
    c_j         = 1.0/((y(jj+1)-y(jj))*(y_cp(jj+1)-y_cp(jj)));
    r(ii,ii)    = -c_j;
end