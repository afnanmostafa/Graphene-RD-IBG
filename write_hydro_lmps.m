function [] = write_hydro_lmps(hydrogenfile,atom_type,x_bot,y_bot,z_bot,z_top,x,y,hydro_C,hydro_C_top,cout,cout2)
%   write list of hydrogens

fid = fopen(hydrogenfile,'w');
l=1;

for m=1:cout
    fprintf(fid,'%g 3 %g %g %g\n',length(atom_type)+l,x_bot(hydro_C(m)), y_bot(hydro_C(m)),z_bot(m)-0.9);
    l=l+1;
end
clear m l
l = 1;
for m=1:cout2
    fprintf(fid,'%g 3 %g %g %g\n',length(atom_type)+cout2+l,x(hydro_C_top(m)), y(hydro_C_top(m)),z_top(m)+0.9);
    l=l+1;
end

fclose(fid);
end

