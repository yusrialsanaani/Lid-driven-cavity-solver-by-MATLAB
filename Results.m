%% Function to plot the results
function Results(W_i,gamma,R,x,y,nx,t,Re)
%% Evaluate primitive variables,[rho,u,v,p,E]
[T,p,u,v]=variables(W_i,gamma,R);

figure(1)
contourf(x,y,sqrt(u'.^2+v'.^2),40);colormap(jet)
xlabel('x','fontweight','bold'); ylabel('y','fontweight','bold');
set(gca,'FontName', 'Times New Roman','FontSize',12,'linewidth',1);
title({['Velocity Field Contour at \itt = ' num2str(t),' and nx = ny = ',num2str(nx)]})

figure(2)
contourf(x,y,p',30); xlabel('x','fontweight','bold'); ylabel('y','fontweight','bold')
set(gca,'FontName', 'Times New Roman','FontSize',12,'linewidth',1);colormap(jet)
title({['Pressure Contour, p at \itt = ' num2str(t),' and nx = ny = ',num2str(nx)]})

figure(3)
quiver(x,y,u',v','r');figure(gcf); axis([0 1 0 1]);
xlabel('x','fontweight','bold');ylabel('y','fontweight','bold');
set(gca,'FontName', 'Times New Roman','FontSize',12,'linewidth',1);
title({['Velocity Vectors at \itt = ' num2str(t),' and nx = ny = ',num2str(nx)]})
%% get the data of the velocity components of Ghia et al.
yref=[0 0.0546875 0.0625 0.0703125 0.1015625 0.171875 0.28125 0.453125 0.5...
    0.6171875 0.734375 0.8515625 0.953125 0.9609375 0.96875 0.9765625 1];
uref=[0 -0.03717 -0.04192 -0.04775 -0.06434 -0.1015 -0.15662 -0.2109...
    -0.20581 -0.13641 0.00332 0.23151 0.68717 0.73722 0.78871 0.84123 1];

xref=[0 0.0625 0.0703 0.0781 0.0938 0.1563 0.2266 0.2344 0.5 0.8047 0.8594 ...
    0.9063 0.9453 0.9531 0.9609 0.9688 1];
vref=[0 0.09233 0.10091 0.10890 0.12317 0.16077 0.17507 0.17527 0.05454 ...
    -0.24533 -0.22445 -0.16914 -0.10313 -0.08864 -0.07391 -0.05906 0];
%% Line Plots for u & v velocities with reference data
figure(4)
u=u';
plot(y,u(:,round(nx/2)),yref,uref,'--','LineWidth',2.0); 
legend('My Result','Ghia Result');legend('boxoff') 
title({['Velocity Component in x at \itt = ' num2str(t),' and nx = ny = ',num2str(nx)]})
xlabel('x','fontweight','bold');ylabel('u','fontweight','bold');
set(gca,'FontName', 'Times New Roman','FontSize',12,'linewidth',1.5,'box','off','XAxisLocation','origin');

figure(5)
plot(x,v(:,round(nx/2)),xref,vref,'--','LineWidth',2.0); 
legend('My Result','Ghia Result');legend('boxoff') 
title({['Velocity Component in y at \itt = ' num2str(t),' and nx = ny = ',num2str(nx)]})
xlabel('x','fontweight','bold'); ylabel('u','fontweight','bold');
set(gca,'FontName', 'Times New Roman','FontSize',12,'linewidth',1.5,'box','off','XAxisLocation','origin');
end
