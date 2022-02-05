%% Function to obtain the total flux (inviscid & Viscus flux)
function [F,G]=Tot_Flux(W,gamma,R,nx,ny,dx,dy,mu,k,Tw)
[T,p,u,v]=variables(W,gamma,R);
[T,p,u,v]=BCs(T,p,u,v,nx,ny,Tw);
[Fi,Gi]=Flux_Inviscid(W,gamma,p); 
[Fv,Gv]=Flux_Viscus(u,v,T,nx,ny,dx,dy,mu,k);
F=Fi+Fv; 
G=Gi+Gv;
end