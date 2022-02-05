%% Function to update Boundary conditions for u,v,p,T
function [T,p,u,v]=BCs(T,p,u,v,nx,ny,Tw)
%% u Boundary conditions
% Down               % Right
u(:,1)=0;            u(nx,:)=0;
% Left               % Up
u(1,:)=0;            u(:,ny)=1;
%% v Boundary conditions
% Down               % Right
v(:,1)=0;            v(nx,:)=0; 
% Left               % Up
v(1,:)=0;            v(:,ny)=0;
%% p Boundary conditions
% Down                         % Right
p(2:nx-1,1)=p(2:nx-1,2);       p(nx,2:ny-1)=p(nx-1,2:ny-1);
% Left                         % Up
p(1,2:ny-1)=p(2,2:ny-1);       p(2:nx-1,ny)=p(2:nx-1,ny-1);
% DL corner          % UL corner
p(1,1)=p(2,2);       p(1,ny)=p(2,ny-1);
% DR corner           % UR corner
p(nx,1)=p(nx-1,2);   p(nx,ny)=p(nx-1,ny-1);
%% T Boundary conditions
% Down               % Right
T(:,1)=Tw;           T(nx,:)=Tw;
% Left               % Up
T(1,:)=Tw;           T(:,ny)=Tw;
end