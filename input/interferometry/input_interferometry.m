%==========================================================================
% input parameters for interferometry
%==========================================================================

%- frequency sampling in Hz -----------------------------------------------
%- The sampling should be evenly spaced for the inverse F transform.-------
%- The sampling must also be sufficiently dense in order to avoid artefacts
%- on the positive time axis in the time-domain source function. This can 
%- be checked with "plot_correlation_source_function".

f_sample=-1000:10:1000;
f_sample=f_sample/5000.0;

%- characteristics of the noise spectrum ----------------------------------

f_peak=1.0/16.0;       % peak frequency in Hz
bandwidth=0.2/16.0;    % bandwidth in Hz
