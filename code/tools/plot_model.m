function plot_model(X,Z,rho,mu)

figure

subplot(1,2,1)
pcolor(X,Z,mu');
axis image
shading flat
title('mu [N/m^2]')
xlabel('x [m]');
ylabel('z [m]');
colorbar
    
subplot(1,2,2)
pcolor(X,Z,rho');
axis image
shading flat
title('rho [kg/m^3]')
xlabel('x [m]');
ylabel('z [m]');
colorbar