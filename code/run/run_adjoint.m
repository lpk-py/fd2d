function K=run_adjoint

%==========================================================================
% run forward simulation
%
% output:
%--------
% K: sensitivity kernel with respect to density
%==========================================================================

%==========================================================================
% read input and set paths
%==========================================================================

input_parameters;
nt=5*round(nt/5);

path(path,'helper_programmes/');

load('cm.mat');
load('v_forward.mat');

%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------

[mu,rho]=define_material_parameters(nx,nz,model_type);

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
    
%- dynamic fields ---------------------------------------------------------

v=zeros(nx,nz);
K=zeros(nx,nz);

sxy=zeros(nx-1,nz);
szy=zeros(nx,nz-1);

%- read adjoint source locations ------------------------------------------

fid=fopen('./seismic_sources/adjoint/source_locations','r');
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

%- read adjoint source time functions -------------------------------------

ns=length(src_x);
stf=zeros(ns,nt);

for n=1:ns
    fid=fopen(['./seismic_sources/adjoint/src_' num2str(n)],'r');
    stf(n,1:nt)=fscanf(fid,'%g',nt);
end

%- compute indices for adjoint source locations ---------------------------

src_x_id=zeros(1,ns);
src_z_id=zeros(1,ns);

x=0:dx:Lx;
z=0:dz:Lz;

for i=1:ns

    src_x_id(i)=min(find(min(abs(x-src_x(i)))==abs(x-src_x(i))));
    src_z_id(i)=min(find(min(abs(z-src_z(i)))==abs(z-src_z(i))));

end
   
%==========================================================================
% iterate
%==========================================================================

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
    
    %- plot and compute kernel every 10th iteration -----------------------
    
    if (mod(n,5)==0)
        
        %- adjoint field --------------------------------------------------
        subplot(2,2,1);
        pcolor(X,Z,v');
        hold on
        for k=1:ns
            plot(src_x(k),src_z(k),'kx');
        end
        hold off
        
        caxis([-0.8*max(max(abs(v))) 0.8*max(max(abs(v)))]);
        colormap(cm);
        axis image
        shading interp
        xlabel('x [m]');
        ylabel('z [m]');
        title('adjoint velocity field [m/s]');
        
        %- forward field --------------------------------------------------
        dummy=squeeze(v_forward(n/5,:,:));
        
        subplot(2,2,2);
        pcolor(X,Z,dummy');
        
        caxis([-0.8*max(max(abs(dummy))) 0.8*max(max(abs(dummy)))]);
        colormap(cm);
        axis image
        shading interp
        xlabel('x [m]');
        ylabel('z [m]');
        title('forward velocity field [m/s]');
        
        %- interaction ----------------------------------------------------
        interaction=v.*dummy;
         
        subplot(2,2,3);
        pcolor(X,Z,interaction');
        
        caxis([-0.8*max(max(abs(interaction))) 0.8*max(max(abs(interaction)))]);
        colormap(cm);
        axis image
        shading interp
        xlabel('x [m]');
        ylabel('z [m]');
        title('interaction (forward \cdot adjoint)');
        
        %- kernel ---------------------------------------------------------
        K=K-interaction*dt; % The minus sign is needed because we run backwards in time.
         
        subplot(2,2,4);
        pcolor(X,Z,K');
        
        caxis([-0.5*max(max(abs(K))) 0.5*max(max(abs(K)))]);
        colormap(cm);
        axis image
        shading interp
        xlabel('x [m]');
        ylabel('z [m]');
        title('sensitivity kernel');
        
        pause(0.01)
        
    end
    
end
