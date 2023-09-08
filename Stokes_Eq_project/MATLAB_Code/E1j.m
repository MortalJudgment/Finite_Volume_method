% Matrix defined gradient on x
function r = E1j(x,Nx,jj)
r = zeros(Nx,Nx);
for ii=1:Nx
    if ii==1
        r(ii,ii) = -1/(2*(x(ii+1)-x(ii)));
        r(ii,ii+1) = 1/(2*(x(ii+1)-x(ii)));
    elseif ii==Nx
        r(ii,ii-1) = -1/(2*(x(ii+1)-x(ii)));
        r(ii,ii) = 1/(2*(x(ii+1)-x(ii)));
    else
        r(ii,ii+1) = 1/(2*(x(ii+1)-x(ii)));
        r(ii,ii-1) = -1/(2*(x(ii+1)-x(ii)));
    end
end

