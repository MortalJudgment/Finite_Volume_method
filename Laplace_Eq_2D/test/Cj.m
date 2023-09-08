function r = Cj(y,y_cp,Nx,jj)
r = zeros(Nx,Nx);
for ii=1:Nx+2
    c_j         = 1.0/((y(jj)-y(jj-1))*(y_cp(jj)-y_cp(jj-1)));
    r(ii,ii)    = -c_j;
end