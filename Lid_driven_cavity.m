%% =========================================================================
%  Matlab program to solve Lid driven cavity.
%  Done by Yusri Al-Sanaani
%% Initializations 
clc; clear all; close all
L=1;        nx=260;     ny=nx;      dx=L/(nx-1);  dy=L/(ny-1);  x=0:dx:L;  
y=0:dy:L;   tmax=20;   t=0;        dt=0.001;     Re=100;       Ma=0.1;  
Pr=0.71;    rhoi=1;    gamma=1.4;  R=1;  U=1;
cp=(gamma*R)/(gamma-1);
Tw=(U/Ma)^2/(gamma*R);
pi=rhoi*R*Tw;
mu=(rhoi*L*U)/Re; 
k=(gamma*mu)/Pr;
sigma=-dt/(2*dx*dy);
%% Obtain initial W=[rho, rho*u, rho*v, rho*E]
u=zeros(nx);        v=zeros(ny); 
p=pi*ones(nx);      rho=rhoi*ones(nx); 
rhoE=(p/(gamma-1))+0.5*rho.*(u.^2+v.^2);
W(:,:,1)=rho;       W(:,:,2)=rho.*u; 
W(:,:,3)=rho.*v;    W(:,:,4)=rhoE;
%% Compute primitive variables         Boundary conditions
[T,p,u,v]=variables(W,gamma,R);        [T,p,u,v]=BCs(T,p,u,v,nx,ny,Tw);
%% Compute the initial flow,  W=[rho, rho*u, rho*v, rhoE]
W=LDC(gamma,u,v,p,T,R);
%% Second-order central FDM scheme in space and  RK4 method in time
W_i=W; 
while t<=tmax
    W1=W_i;
    [F,G]=Tot_Flux(W1,gamma,R,nx,ny,dx,dy,mu,k,Tw);
    R1=RNS(F,G,nx,ny,sigma);
    W2=W_i+(dt/2)*R1;
    [F,G]=Tot_Flux(W2,gamma,R,nx,ny,dx,dy,mu,k,Tw);
    R2=RNS(F,G,nx,ny,sigma);
    W3=W_i+(dt/2)*R2;
    [F,G]=Tot_Flux(W3,gamma,R,nx,ny,dx,dy,mu,k,Tw);
    R3=RNS(F,G,nx,ny,sigma);
    W4=W_i+dt*R3;
    [F,G]=Tot_Flux(W4,gamma,R,nx,ny,dx,dy,mu,k,Tw);
    R4=RNS(F,G,nx,ny,sigma);
    W_RK4=W_i+(dt/6)*(R1+2*R2+2*R3+R4);
    [T,p,u,v]=variables(W_RK4,gamma,R); 
    [T,p,u,v]=BCs(T,p,u,v,nx,ny,Tw);  
    WNEW=LDC(gamma,u,v,p,T,R);
    W_i=WNEW;
    t=t+dt; 
end
%% Reslt Plots
Results(W_i,gamma,R,x,y,nx,t,Re)



