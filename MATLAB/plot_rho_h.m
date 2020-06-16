function [plot_rho,plot_h]=plot_rho_h(rho,thickness)

nl=length(rho);
for m=1:nl
    plot_rho(m*2-1)=rho(m);
    plot_rho(m*2)=rho(m);
end
dmax=10000;
% dmax=maxthick;
% while sum(thickness)>dmax
%     dmax=dmax1000;
% end
plot_h(1)=0;
plot_h(2*nl)=dmax;
h1=0;
for m=1:nl-1
    h1=h1+thickness(m);
    plot_h(m*2)=h1;
    plot_h(m*2+1)=h1;
end

end