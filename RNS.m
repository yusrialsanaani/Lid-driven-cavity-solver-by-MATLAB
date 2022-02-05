%% Function to calculate the RHS of NS equation
function RW=RNS(F,G,nx,ny,sigma)
RF=zeros(size(F));
RG=zeros(size(F));
for n=1:4
    for i=2:nx-1
        for j=2:nx-1
            R11=(F(i+1,j,n)-F(i-1,j,n));
            R22=(G(i,j+1,n)-G(i,j-1,n));
            RF(i,j,n)=R11;
            RG(i,j,n)=R22;
        end
    end
end
    % Three points stencil finite difference schemes
for k=1:4
    for i=1:nx
        RF(1,i,k)=((-3*F(1,i,k)+4*F(2,i,k)-F(3,i,k)));%L
        RF(nx,i,k)=((F(nx-2,i,k)-4*F(nx-1,i,k)+3*F(nx,i,k)));%R
        RG(i,1,k)=((-3*G(i,1,k)+4*G(i,2,k)-G(i,3,k)));%B
        RG(i,ny,k)=((G(i,nx-2,k)-4*G(i,nx-1,k)+3*G(i,ny,k)));%T
    end
end
for n=1:4
    for i=2:nx-1
        RF(i,1,n)=(F(i+1,1,n)-F(i-1,1,n));
        RF(i,ny,n)=(F(i+1,ny,n)-F(i-1,ny,n));
        RG(1,i,n)=(G(1,i+1,n)-G(1,i-1,n));
        RG(nx,i,n)=(G(nx,i+1,n)-G(nx,i-1,n));
    end
end
RW=sigma*(RF+RG);
end
