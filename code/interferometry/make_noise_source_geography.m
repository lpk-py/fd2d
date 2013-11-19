%- compute geographic distribution of noise sources -----------------------

noise_source_distribution=ones(nx,nz);

%noise_source_distribution(120:end,:)=0.0;

%- plot geographic distribution of noise sources --------------------------

figure;
set(gca,'FontSize',20);
load cm_psd

pcolor(X,Z,noise_source_distribution');
shading interp
colormap(cm_psd)
xlabel('x [m]','FontSize',20);
ylabel('z [m]','FontSize',20);
title('power-spectral density distribution of noise sources','FontSize',20);