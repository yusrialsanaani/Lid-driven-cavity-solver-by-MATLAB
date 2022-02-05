%% Function to compute primitive variables
function [T,p,u,v]=variables(W,gamma,R)
u=W(:,:,2)./W(:,:,1);     v=W(:,:,3)./W(:,:,1); 
E=W(:,:,4)./W(:,:,1);     p=(gamma-1)*W(:,:,1).*(E-0.5*(u.^2 +v.^2));
T=p./(R*W(:,:,1));
end