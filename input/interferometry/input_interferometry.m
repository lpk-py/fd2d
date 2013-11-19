%==========================================================================
% input parameters for interferometry
%==========================================================================

%- frequency sampling in Hz -----------------------------------------------
%- The sampling should be evenly spaced for the inverse F transform.-------
%- The sampling must also be sufficiently dense in order to avoid artefacts
%- on the positive time axis in the time-domain source function. This can 
%- be checked with "plot_correlation_source_function".

f_sample=-1000:5:1000;

%- characteristics of the noise spectrum ----------------------------------

f_peak=300.0;       % peak frequency in Hz
bandwidth=100.0;    % bandwidth in Hz