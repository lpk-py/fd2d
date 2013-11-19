%- compute source-time function -------------------------------------------

if strcmp(simulation_mode,'forward_correlation')
    
    stf=1e9*ones(1,length(t));
    
else
    
    stf=zeros(1,length(t));
    stf(1)=1e9;
    stf=butterworth_lp(stf,t,5,f_max,'silent');
    stf=butterworth_hp(stf,t,3,f_min,'silent');
    
end

%- plot -------------------------------------------------------------------

figure
set(gca,'FontSize',20)
hold on

plot(t,stf,'k');
xlabel('time [s]','FontSize',20);
title('source-time function','FontSize',20);
