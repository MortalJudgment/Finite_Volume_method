% Matrix defined gradient on y
function r = E2j(y,Nx,jj)
r = zeros(Nx,Nx);
for ii=1:Nx
    r(ii,ii) = 1/(2*(y(jj+1)-y(jj)));
end