function [tri] = vtk_vector_out(fname, x, y, v, tri)


fid = fopen(fname, 'w');
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Comment\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i float\n', length(x));
for i=1:length(x)
    fprintf(fid,'%g %g 0\n', x(i), y(i));
end
if nargin < 5
    tri = delaunay(x,y)-1;
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
fprintf(fid,'POINT_DATA %i\n', length(x));

fprintf(fid,'SCALARS u_x float 1\n');
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:length(x)
    fprintf(fid,'%g\n', v(i,1));
end
fprintf(fid,'SCALARS u_y float 1\n');
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:length(x)
    fprintf(fid,'%g\n', v(i,2));
end
fprintf(fid,'VECTORS u float\n');
for i=1:length(x)
    fprintf(fid,'%g %g 0\n',v(i,1), v(i,2));
end
fclose(fid);
