%% compute W=[rho, rho*u, rho*v, rhoE]
function W=LDC(gamma,u,v,p,T,R)
rho=p./(R*T);
% E=(p/((gamma-1)*rho))+0.5*(u.^2 + v.^2);
rhoE=(p/(gamma-1))+0.5*rho.*(u.^2 + v.^2);
W(:,:,1)=rho; W(:,:,2)=rho.*u; W(:,:,3)= rho.*v; W(:,:,4)=rhoE;
end