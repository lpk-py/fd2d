if (mod(n,4)==0)
    
    pcolor(X,Z,v')
    axis image
    hold on
    
    if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_correlation'))
        
        for k=1:ns
            plot(src_x(k),src_z(k),'kx')
        end
    end
    
    for k=1:n_receivers
        plot(rec_x(k),rec_z(k),'ko')
    end
    
    hold off
    
    if (n<length(t)/2)
        scale=max(max(abs(v)));
    end
    
    caxis([-scale scale])
    colormap(cm);
    shading interp
    xlabel('x [m]');
    ylabel('z [m]');
    title('velocity field [m/s]');
        
    pause(0.01)
            
end