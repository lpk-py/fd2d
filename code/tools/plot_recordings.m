function plot_recordings(u,t,rec_x,rec_z)

m=max(max(abs(u)));

figure
set(gca,'FontSize',20)
hold on

for k=1:length(rec_x)
    
    plot(t,0.5*m*(k-1)+u(k,:),'k')
    text(min(t),0.5*m*(k-1)+0.2*m,['x=' num2str(rec_x(k)) ', z=' num2str(rec_z(k))],'FontSize',14)
    
end

xlabel('time [s]','FontSize',20);