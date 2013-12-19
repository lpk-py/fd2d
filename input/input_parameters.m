%==========================================================================
% path where seismic sources are located
%==========================================================================

source_path='../../input/sources/forward/';
adjoint_source_path='../../input/sources/adjoint/';

%==========================================================================
% set basic simulation parameters
%==========================================================================

Lx=1e6;     % model extension in x-direction [m]
Lz=1e6;     % model extension in y-direction [m]

nx=200;     % grid points in x-direction
nz=200;     % grid points in z-direction

dt=0.5;     % time step [s]
nt=400;     % number of iterations

order=4;    % finite-difference order (2 or 4)

%==========================================================================
% model type
%==========================================================================

model_type=1;

% 1=homogeneous 
% 2=homogeneous with localised density perturbation
% 3=layered medium
% 4=layered with localised density perturbation
% 5=vertical gradient medium
% 6=vertical gradient medium with localised density perturbation

% 'initial'= read initial model for waveform inversion (mu_initial, rho_initial)

%==========================================================================
% source-time function
%==========================================================================

f_min=0.02;     % minimum frequency [Hz]
f_max=0.10;     % maximum frequency [Hz]

%==========================================================================
% simulation mode
%==========================================================================

simulation_mode='forward_correlation';

% 'forward'                 regular forward simulation
% 'forward_correlation'     forward simulation where Fourier transform is computed on-the-fly (preparation to compute correlation function)
% 'correlation'             compute correlation functions as output, requires Fourier-transformed Green function to be present
% 'noise_source_kernel'     compute sensitivity kernel for the noise source power-spectral density distribution

%==========================================================================
% receiver positions
%==========================================================================

%- a large number of receivers in an open rectangular configuration
%rec_x=[50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0  50.0  50.0  50.0  60.0 70.0 80.0 90.0 100.0 110.0 120.0 130.0 60.0  70.0  80.0  90.0  100.0 110.0 120.0 130.0];
%rec_z=[70.0 80.0 90.0 100.0 110.0 120.0 130.0 140.0 150.0 160.0 170.0 180.0 70.0 70.0 70.0 70.0 70.0  70.0  70.0  70.0  180.0 180.0 180.0 180.0 180.0 180.0 180.0 180.0];

%- just one receiver
rec_x=[750000.0];
rec_z=[500000.0];

%- a large number of receivers in a closed rectangular configuration
%rec_x=[50.0  50.0  50.0  50.0  50.0   50.0    70.0  90.0 110.0 130.0   70.0  90.0 110.0 130.0  150.0 150.0 150.0 150.0 150.0  150.0];
%rec_z=[70.0  90.0 110.0 130.0 150.0  170.0    70.0  70.0  70.0  70.0  170.0 170.0 170.0 170.0   70.0  90.0 110.0 130.0 150.0  170.0];

%rec_x=zeros(1,49);
%rec_z=zeros(1,49);
%n=1;

%for i=200e3:100e3:800e3
%    for j=200e3:100e3:800e3
%        rec_x(n)=i;
%        rec_z(n)=j;
%        n=n+1;
%    end
%end

%==========================================================================
% absorbing boundaries
%==========================================================================

width=100000.0;     % width of the boundary layer in km

absorb_left=1;  % absorb waves on the left boundary
absorb_right=1; % absorb waves on the right boundary
absorb_top=1;   % absorb waves on the top boundary
absorb_bottom=1;% absorb waves on the bottom boundary

%==========================================================================
% make wavepropagation movie
%==========================================================================

make_movie='no';                           % 'yes' or 'no'
movie_file='../../output/cc_onesided_distribution.mp4';    % output file name, should be .mp4