function misfit=make_adjoint_sources(u,u_pert,t)

%==========================================================================
% compute and store adjoint sources
%
% input:
%-------
% u: recordings for the original (unperturbed) medium
% u_pert: recordings for the perturbed medium
% t: time axis
% rec_x, rec_z: receiver positions in x- and z-directions
%==========================================================================

%==========================================================================
%- initialisations --------------------------------------------------------
%==========================================================================

path(path,'../input/');
input_parameters;

fid_loc=fopen([adjoint_source_path 'source_locations'],'w');

nt=length(t);

misfit=0.0;

%==========================================================================
%- march through the various recodings ------------------------------------
%==========================================================================

for n=1:length(rec_x)
   
    fprintf(1,'station number %d\n',n)
    
    %- plot traces --------------------------------------------------------
    
    plot(t,u(n,:),'k')
    hold on
    plot(t,u_pert(n,:),'r')
    plot(t,u(n,:)-u_pert(n,:),'k--')
    hold off
   
    title(['receiver ' num2str(n) ' ,original in black, perturbed in red, difference dashed'])
    xlabel('t [s]')
    ylabel('displacement [m]')
   
    %- select time windows ------------------------------------------------
   
    disp('select left window');
    [left,dummy]=ginput(1);
    disp('select_right_window');
    [right,dummy]=ginput(1);
   
    %- taper seismogram difference and compute misfit ---------------------
   
    adstf=taper(u(n,:)-u_pert(n,:),t,left,right,(right-left)/10);
    misfit=misfit+sum(adstf.*adstf)*(t(2)-t(1));
   
    plot(t,adstf,'k')
    xlabel('t [s]')
    title('adjoint source before time reversal')
    pause(1.0)
   
    %- write time-reversed adjoint source to file -------------------------
   
    fprintf(fid_loc,'%g %g\n',rec_x(n),rec_z(n));
    
    %- write source time functions ------------------------------------
    fn=[adjoint_source_path 'src_' num2str(n)];
    fid_src=fopen(fn,'w');
    for k=1:nt
        fprintf(fid_src,'%g\n',adstf(nt+1-k));
    end
    fclose(fid_src);
      
end

%==========================================================================
%- clean up ---------------------------------------------------------------
%==========================================================================

fclose(fid_loc);