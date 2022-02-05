%% Function to Evaluate the invicid flux
function [Fi,Gi]=Flux_Inviscid(W,gamma,p)
%% Conserved variables
rho=W(:,:,1); rhou=W(:,:,2); rhov=W(:,:,3); E=W(:,:,4)./rho;
%% compute primitive variables
u=rhou./rho; v=rhov./rho; 
% p=(gamma-1)*rho.*(E-0.5*(u.^2 +v.^2));
%% Compute flux functions
Fi=zeros(size(W));
Fi(:,:,1)=rhou; Fi(:,:,2)=rhou.*u+p; Fi(:,:,3)=rhov.*u; Fi(:,:,4)=u.*(rho.*E+p);
Gi=zeros(size(W));
Gi(:,:,1)=rhov; Gi(:,:,2)=rhou.*v; Gi(:,:,3)=rhov.*v+p; Gi(:,:,4)=v.*(rho.*E+p);
end