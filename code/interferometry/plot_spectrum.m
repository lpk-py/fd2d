function plot_spectrum(x_plot,z_plot)

%- set paths --------------------------------------------------------------

path(path,'../../input/');
path(path,'../../input/interferometry');

%- load the displacement spectrum field -----------------------------------

load('../../output/G.mat');

%- read input and make coordinates ----------------------------------------

input_parameters;
input_interferometry;

dx=Lx/(nx-1);
dz=Lz/(nz-1);

x=0:dx:Lx;
z=0:dz:Lz;

%- find index -------------------------------------------------------------

x_id=min(find(min(abs(x-x_plot))==abs(x-x_plot)));
z_id=min(find(min(abs(z-z_plot))==abs(z-z_plot)));

%- get the spectrum -------------------------------------------------------

s=reshape(G(x_id,z_id,:),length(f_sample),1);

%- plot results -----------------------------------------------------------

figure

subplot(2,1,1)
plot(f_sample,abs(s),'k');
xlabel('\nu [Hz]')
ylabel('amplitude spectrum')

title(['velocity spectrum at position x=' num2str(x(x_id)) ' m, z=' num2str(z(z_id)) ' m'])

subplot(2,1,2)
plot(f_sample,angle(s),'k');
xlabel('\nu [Hz]')
ylabel('phase spectrum')

