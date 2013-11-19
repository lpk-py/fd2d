
%- compute noise spectrum -------------------------------------------------

noise_spectrum=exp(-(abs(f_sample)-f_peak).^2/bandwidth^2);

%- plot -------------------------------------------------------------------

figure
set(gca,'FontSize',20)
hold on

plot(f_sample,noise_spectrum,'k');
xlabel('frequency [Hz]','FontSize',20);
title('noise power-spectral density','FontSize',20);
