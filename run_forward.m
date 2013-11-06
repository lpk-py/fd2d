function [u,t,rec_x,rec_z]=run_forward(plot_mode)

%==========================================================================
% run forward simulation
%
% input:
%-------
% plot_mode: show wavefield simulation when plot_mode=='plot'
%
% output:
%--------
% u: displacement seismograms
% t: time axis
% rec_x: x-coordinate of receiver positions
% rec_z: z-coordinate of receiver positions
%
%==========================================================================

%==========================================================================
% read input and set paths
%==========================================================================

input_parameters;

path(path,'helper_programmes/');

if strcmp(plot_mode,'plot')
    h_vel=figure;
    load cm;
end

%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------

[mu,rho]=define_material_parameters(nx,nz);

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

if strcmp(plot_mode,'plot')
    
    h_model=figure;
    
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
        
end
    
%- dynamic fields ---------------------------------------------------------

v=zeros(nx,nz);

sxy=zeros(nx-1,nz);
szy=zeros(nx,nz-1);

%- recommended time step from CFL criterion -------------------------------

%vmax=max(max(sqrt(mu./rho)));
%h=min(dx,dz);

%fprintf(1,'recommended time step: dt=%g s\n',0.5*h/vmax);
%fprintf(1,'used time step: dt=%g s\n',dt);

%- read seismic source locations ------------------------------------------

fid=fopen([source_path 'source_locations'],'r');
src_x=zeros(1);
src_z=zeros(1);

k=1;
while (feof(fid)==0)
    src_x(k)=fscanf(fid,'%g',1);
    src_z(k)=fscanf(fid,'%g',1);
    fgetl(fid);
    k=k+1;
end

fclose(fid);

%- read source time functions ---------------------------------------------

ns=length(src_x);
stf=zeros(ns,nt);

for n=1:ns
    fid=fopen([source_path 'src_' num2str(n)],'r');
    stf(n,1:nt)=fscanf(fid,'%g',nt);
end

%- compute indices for source locations -----------------------------------

src_x_id=zeros(1,ns);
src_z_id=zeros(1,ns);

x=0:dx:Lx;
z=0:dz:Lz;

for i=1:ns

    src_x_id(i)=min(find(min(abs(x-src_x(i)))==abs(x-src_x(i))));
    src_z_id(i)=min(find(min(abs(z-src_z(i)))==abs(z-src_z(i))));

end

%- compute indices for receiver locations ---------------------------------

n_receivers=length(rec_x);

rec_x_id=zeros(1,n_receivers);
rec_z_id=zeros(1,n_receivers);

x=0:dx:Lx;
z=0:dz:Lz;

for i=1:n_receivers

    rec_x_id(i)=min(find(min(abs(x-rec_x(i)))==abs(x-rec_x(i))));
    rec_z_id(i)=min(find(min(abs(z-rec_z(i)))==abs(z-rec_z(i))));

end
   
%==========================================================================
% initialise seismograms
%==========================================================================

displacement_seismograms=zeros(n_receivers,nt);
velocity_seismograms=zeros(n_receivers,nt);

%==========================================================================
% iterate
%==========================================================================

t=0;

figure(h_vel);

for n=1:nt
    
    %- compute divergence of current stress tensor and add external forces
    
    DS=div_s(sxy,szy,dx,dz,nx,nz,order);
    
    for i=1:ns
        DS(src_x_id(i),src_z_id(i))=DS(src_x_id(i),src_z_id(i))+stf(i,n);
    end
    
    %- update velocity field ----------------------------------------------
    
    v=v+dt*DS./rho;
    
    %- compute derivatives of current velocity and update stress tensor ---
    
    sxy=sxy+dt*mu(1:nx-1,1:nz).*dx_v(v,dx,dz,nx,nz,order);
    szy=szy+dt*mu(:,1:nz-1).*dz_v(v,dx,dz,nx,nz,order);
    
    %- record velocity seismograms ----------------------------------------
    
    for k=1:n_receivers
        velocity_seismograms(k,n)=v(rec_x_id(k),rec_z_id(k));
    end
    
    %- update time --------------------------------------------------------
    
    t=t+dt;
    
    %- plot velocity field ------------------------------------------------
    
    if strcmp(plot_mode,'plot')
    
        pcolor(X,Z,v');
    
        m=max(max(abs(v)));
        caxis([-m,m]);
        colormap(cm);
        axis image
        shading interp
        xlabel('x [m]');
        ylabel('z [m]');
        title('velocity field [m/s]');
        
        hold on
        for k=1:ns
            plot(src_x(k),src_z(k),'kx')
        end
        for k=1:n_receivers
            plot(rec_x(k),rec_z(k),'ko')
        end 
        hold off

        pause(0.02)
        
    end
    
end

%==========================================================================
% output
%==========================================================================

t=0:dt:dt*(nt-1);
u=cumsum(velocity_seismograms,2);
u(:,nt)=u(:,nt-1);