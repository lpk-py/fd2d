function make_sources(frequency, src_x, src_z)

%==========================================================================
% write seismic sources with Ricker wavelet as source-time function
%
% input:
%-------
% frequency: dominant frequency of the Ricker wavelet
% src_x, src_z: locations of the sources in x- and z-directions
%==========================================================================

%- initialisations --------------------------------------------------------

path(path,'helper_programmes/');
input_parameters;

nt=5*round(nt/5);

n_src=length(src_x);

t=0:dt:dt*(nt-1);

%- write sources ----------------------------------------------------------

fid_loc=fopen('seismic_sources/forward/source_locations','w');

for n=1:n_src
    
    %- write source position ----------------------------------------------
    fprintf(fid_loc,'%g %g\n',src_x(n),src_z(n));
    
    %- write source time functions ----------------------------------------
    fn=['seismic_sources/forward/src_' num2str(n)];
    fid_src=fopen(fn,'w');
    s=1e12*ricker_wavelet(1/(2*pi*frequency),1/(1.5*frequency),t);
    %- The factor 1e12 ensures that the wave field is not too small.
    for k=1:nt
        fprintf(fid_src,'%g\n',s(k));
    end
    
end

fclose(fid_loc);