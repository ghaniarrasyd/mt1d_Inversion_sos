function [misfit]=misfitMT(rho_obs,phase_obs,rho_cal,phase_cal)

d2r=pi/180;
nd=length(rho_obs);
for k=1:nd
    m(k)=(abs(log10(rho_cal(k)/rho_obs(k)))+abs(d2r*phase_cal(k)-d2r*phase_obs(k)));
end
misfit=sum(m)/(nd);
end