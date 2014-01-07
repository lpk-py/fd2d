%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generalities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- number of noise sources
n_noise_sources=2;

%- characteristics of the noise spectrum ----------------------------------
%- only needed in this routine --------------------------------------------
f_peak=1.0/16.0;       % peak frequency in Hz
bandwidth=0.35/16.0;    % bandwidth in Hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% source spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise_spectrum=zeros(length(f_sample),n_noise_sources);

for i=1:n_noise_sources
    
    %- compute noise spectrum for source i --------------------------------
    noise_spectrum(:,i)=0.7*exp(-(abs(f_sample)-0.07).^2/bandwidth^2);
    noise_spectrum(:,i)=noise_spectrum(:,i)+exp(-(abs(f_sample)-0.12).^2/(0.5*bandwidth)^2)';
    
    %- plot noise spectrum for source i -----------------------------------
    figure
    set(gca,'FontSize',20)
    hold on

    plot(f_sample,noise_spectrum(:,i),'k');
    xlabel('frequency [Hz]','FontSize',20);
    title(['noise power-spectral density for source region ', num2str(i)] ,'FontSize',20);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geographic power-spectral density distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise_source_distribution=zeros(nx,nz,n_noise_sources);

for i=1:n_noise_sources
    
    %- compute noise source distribution for each source region -----------
    noise_source_distribution(:,:,i)=1.0;
    
    %- plot noise source distribution -------------------------------------

    figure;
    set(gca,'FontSize',20);
    load cm_psd

    pcolor(X,Z,noise_source_distribution(:,:,i)');
    shading interp
    colormap(cm_psd)
    xlabel('x [m]','FontSize',20);
    ylabel('z [m]','FontSize',20);
    title(['power-spectral density distribution of noise sources for region ', num2str(i)],'FontSize',20);
    
end