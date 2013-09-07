function vtk_output(filename, folder, step, endstep)

points = load(sprintf('%s_points',filename));

for i=0:step:endstep
    fnameout = sprintf('%s/field%i.vtk',folder,i);
    fnamein = sprintf('%s_fieldu_%i.field',filename,i);
    uh_rec = load(fnamein);
    if i==0
	tri = vtk_vector_out(fnameout,points(:,1), points(:,2), uh_rec);
    else
	vtk_vector_out(fnameout,points(:,1), points(:,2), uh_rec, tri);
    end

end


