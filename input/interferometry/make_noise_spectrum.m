
%- compute noise spectrum -------------------------------------------------

%noise_spectrum=exp(-(abs(f_sample)-f_peak).^2/bandwidth^2);

noise_spectrum=0.7*exp(-(abs(f_sample)-0.07).^2/bandwidth^2);
noise_spectrum=noise_spectrum+exp(-(abs(f_sample)-0.14).^2/bandwidth^2);

%- plot -------------------------------------------------------------------

figure
set(gca,'FontSize',20)
hold on

plot(f_sample,noise_spectrum,'k');
xlabel('frequency [Hz]','FontSize',20);
title('noise power-spectral density','FontSize',20);
