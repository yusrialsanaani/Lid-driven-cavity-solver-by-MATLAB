%% Function to obtain the viscus flux
function [Fv,Gv]=Flux_Viscus(u,v,T,nx,ny,dx,dy,mu,k)

%% Find du/dx,du/dy,dv/dx,dv/dy
dudx=zeros(nx);  dudy=zeros(nx);
dvdy=zeros(nx);  dvdx=zeros(nx);
for i=2:nx-1
    for j=2:ny-1
        dudx(i,j)=(u(i+1,j)-u(i-1,j))/(2*dx);
        
        dudy(i,j)=(u(i,j+1)-u(i,j-1))/(2*dy);
        
        dvdx(i,j)=(v(i+1,j)-v(i-1,j))/(2*dx);
        
        dvdy(i,j)=(v(i,j+1)-v(i,j-1))/(2*dy);
    end
end
%%
for i=2:nx-1
    dudx(i,1)=(u(i+1,1)-u(i-1,1))/(2*dx);  % Down
    dudx(i,ny)=(u(i+1,ny)-u(i-1,ny))/(2*dx);  % Up
    dudy(1,i)=(u(1,i+1)-u(1,i-1))/(2*dy);% L
    dudy(nx,i)=(u(nx,i+1)-u(nx,i-1))/(2*dy);% T
    
    dvdx(i,1)=(v(i+1,1)-v(i-1,1))/(2*dx);  % Down
    dvdx(i,ny)=(v(i+1,ny)-v(i-1,ny))/(2*dx);  % Up
    dvdy(1,i)=(v(1,i+1)-v(1,i-1))/(2*dy);% L
    dvdy(nx,i)=(v(nx,i+1)-v(nx,i-1))/(2*dy);% T
end
% Three points one-sided stencil finite difference schemes at walls
for i=1:nx
    dudx(1,i)=(-3*u(1,i)+4*u(2,i)-u(3,i))/(2*dx);         % Left
    dudx(nx,i)=(u(nx-2,i)-4*u(nx-1,i)+3*u(nx,i))/(2*dx);  % Right 
    dudy(i,1)=(-3*u(i,1)+4*u(i,2)-u(i,3))/(2*dx);         % b
    dudy(i,ny)=(u(i,ny-2)-4*u(i,ny-1)+3*u(i,ny))/(2*dx);  % T
    
    dvdx(1,i)=(-3*v(1,i)+4*v(2,i)-v(3,i))/(2*dx);         % Left
    dvdx(nx,i)=(v(nx-2,i)-4*v(nx-1,i)+3*v(nx,i))/(2*dx);  % Right
    dvdy(i,1)=(-3*v(i,1)+4*v(i,2)-v(i,3))/(2*dx);         % B
    dvdy(i,ny)=(v(i,ny-2)-4*v(i,ny-1)+3*v(i,ny))/(2*dx);  % T
end
%% find stresses
ss_xx=(2*mu/3)*(2*dudx-dvdy);
ss_yy=(2*mu/3)*(2*dvdy-dudx);
ss_yx=mu*(dudy+dvdx);
ss_xy=ss_yx;
%% Find Qx=-k*dT/dx and Qy=-k*dT/dy
Qx=zeros(nx);
Qy=zeros(nx);
for i=2:nx-1
    for j=2:ny-1
        Qx(i,j)=(-k/(2*dx))*(T(i+1,j)-T(i-1,j));
        Qy(i,j)=(-k/(2*dy))*(T(i,j+1)-T(i,j-1));
    end
end
%%
for i=2:nx-1
    Qx(i,1)=(-k/(2*dx))*(T(i+1,1)-T(i-1,1));  % Down
    Qx(i,ny)=(-k/(2*dx))*(T(i+1,ny)-T(i-1,ny));  % Up
    Qy(1,i)=(-k/(2*dx))*(T(1,i+1)-T(1,i-1));% L
    Qy(nx,i)=(-k/(2*dx))*(T(nx,i+1)-T(nx,i-1));% T
end
%%

% Three points stencil finite difference schemes
for i=1:nx
    Qx(1,i)=(-k/(2*dx))*(-3*T(1,i)+4*T(2,i)-T(3,i));
    Qx(nx,i)=(-k/(2*dx))*(T(nx-2,i)-4*T(nx-1,i)+3*T(nx,i));
    Qy(i,1)=(-k/(2*dy))*(-3*T(i,1)+4*T(i,2)-T(i,3));
    Qy(i,ny)=(-k/(2*dy))*(T(i,nx-2)-4*T(i,nx-1)+3*T(i,ny));
end
%% Calculate viscus flux
Fv=zeros(nx,ny,4);
Fv(:,:,1)=zeros(nx,ny); Fv(:,:,2)=-ss_xx; Fv(:,:,3)=-ss_xy; Fv(:,:,4)=Qx-u.*ss_xx-v.*ss_yx;
Gv=zeros(nx,ny,4);
Gv(:,:,1)=zeros(nx,ny); Gv(:,:,2)=-ss_xy; Gv(:,:,3)=-ss_yy; Gv(:,:,4)=Qy-u.*ss_xy-v.*ss_yy;
end