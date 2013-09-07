function  vtk_mesh_out(femregion, dir, varargin)
%usage: vtk_mesh_out(filename, coords, t, 'param_name', param , 'param_name', param, ...)

fname = [dir,'/mesh_and_params.vtk'];
x     = femregion.coord(:,1);
y     = femregion.coord(:,2);
tri   = femregion.connectivity(1:3,:)'-1;

fid = fopen(fname, 'w');

fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Comment\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i float\n', length(x));
for i=1:length(x)
    fprintf(fid,'%g %g 0\n', x(i), y(i));
end
numtria = length(tri(:,1));
fprintf(fid,'CELLS %i %i\n',numtria, 4*numtria);
for i=1:numtria
    fprintf(fid,'3 %i %i %i\n', tri(i,1), tri(i,2), tri(i,3));
end
fprintf(fid,'CELL_TYPES %i\n', numtria);
for i=1:numtria
    fprintf(fid,'5\n' );
end

fprintf(fid,'CELL_DATA %i\n', numtria);
for i = 1:2:nargin-3
    par  = varargin{i+1};
    fprintf(fid,'SCALARS %s float 1\n', varargin{i});
    fprintf(fid,'LOOKUP_TABLE default\n');
    if length(par)>1
	for j=1:length(par)
	    fprintf(fid,'%g\n', par(j));
	end
    else
	for j=1:numtria
	    fprintf(fid,'%g\n',par);
	end
    end
end
fclose(fid);
