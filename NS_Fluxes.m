%% Evaluate primitive variables and NS flux functions,[rho,u,v,p,E,F,G]
function [F,G,rho,u,v,p,E]=NS_Fluxes(W,gamma,nx,ny,dx,dy,vis,Tw,k)
% Conserved variables
rho=W(:,:,1); rhou=W(:,:,2); rhov=W(:,:,3); E=W(:,:,4)./rho;
% compute primitive variables
u=rhou./rho; v=rhov./rho; p=(gamma-1)*rho.*(E-0.5*(u.^2 +v.^2));
% u(:,1)=0; u(:,nx)=0; u(1,:)=0; u(ny,:)=1;
% v(:,1)=0; v(:,nx)=0; v(1,:)=0; v(ny,:)=0;
% compute flux functions
%% Find du/dx,du/dy,dv/dx,dv/dy
dudx=zeros(nx,nx);
dudy=zeros(nx,nx);
dvdy=zeros(nx,nx);
dvdx=zeros(nx,nx);
for i=2:nx-1
    for j=2:ny-1
        dudx(j,i)=(u(i+1,j)-u(i-1,j))/(2*dx);
        dudy(j,i)=(u(i,j+1)-u(i,j-1))/(2*dy);
        dvdx(j,i)=(v(i+1,j)-v(i-1,j))/(2*dx);
        dvdy(j,i)=(v(i,j+1)-v(i,j-1))/(2*dy);
    end
end
% Three points stencil finite difference schemes
for i=1:nx
    dudx(1,i)=(-3*u(1,i)+4*u(2,i)-u(3,i))/(2*dx);
    dudx(nx,i)=(u(nx-2,i)-4*u(nx-1,i)+3*u(nx,i))/(2*dx);
    dudy(i,1)=(-3*u(i,1)+4*u(i,2)-u(i,3))/(2*dy);
    dudy(i,ny)=(u(i,nx-2)-4*u(i,nx-1)+3*u(i,ny))/(2*dy);
end
%% find stresses
ss_xx=(2*vis/3)*(2*dudx-dvdy);
ss_yy=(2*vis/3)*(2*dvdy-dudx);
ss_yx=vis*(dudy+dvdx);
ss_xy=ss_yx;
% Find Qx=-k*dT/dx and Qy=-k*dT/dy
T=zeros(ny,nx);
T(:,1)=Tw; T(:,nx)=Tw; T(1,:)=Tw; T(ny,:)=Tw;
Qx=zeros(ny,nx);
Qy=zeros(ny,nx);
for i=2:nx-1
    for j=2:ny-1
        Qx(i,j)=(-k/(2*dx))*(T(i+1,j)-T(i-1,j));
        Qy(i,j)=(-k/(2*dy))*(T(i,j+1)-T(i,j-1));
    end
end
% Three points stencil finite difference schemes
for i=1:nx
    Qx(1,i)=(-k/(2*dx))*(-3*T(1,i)+4*T(2,i)-T(3,i));
    Qx(nx,i)=(-k/(2*dx))*(T(nx-2,i)-4*T(nx-1,i)+3*T(nx,i));
    Qy(i,1)=(-k/(2*dy))*(-3*T(i,1)+4*T(i,2)-T(i,3));
    Qy(i,ny)=(-k/(2*dy))*(T(i,nx-2)-4*T(i,nx-1)+3*T(i,ny));
end

Fi=zeros(size(W)); 
Fi(:,:,1)=rhou; Fi(:,:,2)=rhou.*u+p; Fi(:,:,3)=rhov.*u; Fi(:,:,4)=u.*(rho.*E+p);
Fv=zeros(size(W));
Fv(:,:,1)=zeros(nx,ny); Fv(:,:,2)=-ss_xx; Fv(:,:,3)=-ss_xy; Fv(:,:,4)=Qx-u.*ss_xx-v.*ss_yx;
F=Fi+Fv;
Gi=zeros(size(W));
Gi(:,:,1)=rhov; Gi(:,:,2)=rhou.*v; Gi(:,:,3)=rhov.*v+p; Gi(:,:,4)=v.*(rho.*E+p);
Gv=zeros(size(W));
Gv(:,:,1)=zeros(nx,ny); Gv(:,:,2)=-ss_xy; Gv(:,:,3)=-ss_yy; Gv(:,:,4)=Qy-u.*ss_xy-v.*ss_yy;
G=Gi+Gv;
end