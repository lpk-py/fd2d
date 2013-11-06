%- write time-reversed source positions

fid=fopen('./seismic_sources/reverse/source_locations','w');

for k=1:length(rec_x)
    fprintf(fid,'%g %g\n',rec_x(k),rec_z(k));
end 

fclose(fid);

%- write time-reversed recordings

nt=length(u(1,:));

for k=1:length(rec_x)
    
    filename=['./seismic_sources/reverse/src_' num2str(k)];
    fid=fopen(filename,'w');
    
    for i=0:nt-1
        fprintf(fid,'%g\n',u(k,nt-i));
    end
    
    fclose(fid);
    
end
