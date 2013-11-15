function plot_correlation_source_function(x_plot,z_plot)

%- set paths --------------------------------------------------------------

path(path,'../../input/');
path(path,'../../input/interferometry');

%- load the displacement spectrum field -----------------------------------

load('../../output/G.mat');

%- read input and make space-time coordinates -----------------------------

input_parameters;
input_interferometry;

dx=Lx/(nx-1);
dz=Lz/(nz-1);

x=0:dx:Lx;
z=0:dz:Lz;

t=-(nt-1)*dt:dt:(nt-1)*dt;

%- find index -------------------------------------------------------------

x_id=min(find(min(abs(x-x_plot))==abs(x-x_plot)));
z_id=min(find(min(abs(z-z_plot))==abs(z-z_plot)));

%- get the spectrum -------------------------------------------------------

s=reshape(G(x_id,z_id,:),length(f_sample),1);
s=conj(s);

%- approximate inverse Fourier transform ----------------------------------

stf=zeros(1,length(t));
w_sample=2*pi*f_sample;
dw=w_sample(2)-w_sample(1);

for k=1:length(f_sample)
    stf=stf+s(k)*exp(sqrt(-1)*w_sample(k)*t);
end

stf=dw*stf/(2*pi);

%- plot results -----------------------------------------------------------

%- frequency domain

figure
set(gca,'FontSize',20)

subplot(2,1,1)
plot(f_sample,abs(s),'k');
xlabel('\nu [Hz]','FontSize',20)
ylabel('amplitude spectrum','FontSize',20)

title(['velocity spectrum at position x=' num2str(x(x_id)) ' m, z=' num2str(z(z_id)) ' m'],'FontSize',20)

subplot(2,1,2)
plot(f_sample,angle(s),'k');
xlabel('\nu [Hz]','FontSize',20)
ylabel('phase spectrum','FontSize',20)

%- time domain

figure
set(gca,'FontSize',20)
hold on

plot(t,real(stf),'k')
plot(t,imag(stf),'k--')
xlabel('time [s]','FontSize',20)
ylabel('amplitude','FontSize',20)

title('time-domain source function (solid: real, dashed: imaginary)','FontSize',20)