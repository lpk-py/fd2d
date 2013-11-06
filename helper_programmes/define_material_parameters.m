function [mu,rho]=define_material_parameters(nx,nz)

%==========================================================================
% generate material parameters mu [N/m^2] and rho [kg/m^3]
%
% input: grid points of the velocity and density field in x-direction
% (nx) and z-direction (nz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1=homogeneous, 2=homogeneous with pointwise perturbation

model=2;

%==========================================================================

if (model==1)
    
    rho=3000.0*ones(nx,nz);
    mu=4.8e10*ones(nx,nz);
    
elseif (model==2)
    
    rho=3000.0*ones(nx,nz);
    mu=4.8e10*ones(nx,nz);
    
    rho(100,75)=1000.0;
    
end
    

