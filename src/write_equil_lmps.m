function [] = write_equil_lmps(lmp_input,rand_bot_atoms,rand_top_atoms,len,wid,stacking,bond_density)
%   writes LAMMPS input script for equilibration

lmp_file = sprintf('gr%dx%d%s_%.2fBD.data',len,wid,stacking,bond_density);
fid2 = fopen(lmp_input, 'w');

fprintf(fid2, '#Interlayer bonded Bilayer Graphene with random isolated bonds + passivation by H');
fprintf(fid2,'\n');
fprintf(fid2,'\n');
fprintf(fid2,'dimension 3');
fprintf(fid2,'\n');
fprintf(fid2,'units metal');
fprintf(fid2,'\n');
fprintf(fid2,'processors * 1 1');
fprintf(fid2,'\n');
fprintf(fid2,'boundary p p p');
fprintf(fid2,'\n');
fprintf(fid2,'neighbor 0.3 bin');
fprintf(fid2,'\n');
fprintf(fid2,'neigh_modify every 1 delay 0 check yes');
fprintf(fid2,'\n\n');

fprintf(fid2,'# ======= define variables ======= #');
fprintf(fid2,'\n');
fprintf(fid2,'\n');
fprintf(fid2,'variable T equal 300.0');
fprintf(fid2,'\n');
fprintf(fid2,'variable V equal lx*ly*6.7  # 3.35 (1st layer) + gap (0) + 3.35 (2nd layer)');
fprintf(fid2,'\n');
fprintf(fid2,'variable dt equal 0.0001');
fprintf(fid2,'\n\n');

fprintf(fid2,'# ======= read data ======= #');
fprintf(fid2,'\n');
fprintf(fid2,'\n');
fprintf(fid2,'atom_style atomic');
fprintf(fid2,'\n');
fprintf(fid2,'read_data %s',lmp_file);
fprintf(fid2,'\n\n');

fprintf(fid2,'# ======= potential information ======= #');
fprintf(fid2,'\n');
fprintf(fid2,'\n');
fprintf(fid2,'pair_style			airebo 3 1 0'); %LJ on (1), torsion off (0)
fprintf(fid2,'\n');
fprintf(fid2,'pair_coeff 			* * CH.airebo C C H');
fprintf(fid2,'\n');
fprintf(fid2,'dump 1 all xyz 1000 bondnew.xyz');
fprintf(fid2,'\n');
fprintf(fid2,'dump_modify 1 sort id');
fprintf(fid2,'\n\n');

fprintf(fid2,'# ======= time step ======= #');
fprintf(fid2,'\n');
fprintf(fid2,'\n');
fprintf(fid2,'timestep ${dt}');
fprintf(fid2,'\n\n');

fprintf(fid2,'# ======= bond formation ======= #');
fprintf(fid2,'\n');
fprintf(fid2,'\n');
fprintf(fid2,'group atom1 id \t');
fprintf(fid2, '%d ',rand_bot_atoms);
fprintf(fid2,'\n\n');
fprintf(fid2,'group atom2 id \t');
fprintf(fid2, '%d ',rand_top_atoms);
fprintf(fid2,'\n\n');
fprintf(fid2,'displace_atoms atom1 move 0 0 1.15 units box');
fprintf(fid2,'\n');
fprintf(fid2,'displace_atoms atom2 move 0 0 -1.15 units box');
fprintf(fid2,'\n\n');
fprintf(fid2,'# ======= grouping bonded and non-bonded atoms ======= #\n\n');
fprintf(fid2,'group bonded union atom1 atom2');
fprintf(fid2,'\n');
fprintf(fid2,'group nonbonded subtract all bonded');
fprintf(fid2,'\n');
fprintf(fid2,'\n');

fprintf(fid2,'# ======= thermodynamics output setting ======= #\n\n');
fprintf(fid2,'thermo_style    	custom step time temp press pxx pyy pzz pe ke etotal');
fprintf(fid2,'\n');
fprintf(fid2,'thermo          	1000');
fprintf(fid2,'\n\n');

fprintf(fid2,'# ======= minimization ======= #');
fprintf(fid2,'\n');
fprintf(fid2,'\n');
fprintf(fid2,'min_style		cg');
fprintf(fid2,'\n');
fprintf(fid2,'fix 			1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz');
fprintf(fid2,'\n');
fprintf(fid2,'minimize		0 1e-10 200000 400000\n');
fprintf(fid2,'unfix			1');
fprintf(fid2,'\n\n');

fprintf(fid2,'# ======= equilibration ======= #\n\n');

fprintf(fid2,'fix 2 all recenter INIT INIT INIT');
fprintf(fid2,'\n\n');
fprintf(fid2,'fix 1 nonbonded npt temp $T $T 0.1 x 0.0 0.0 1.0 y 0.0 0.0 1.0 couple xy');
fprintf(fid2,'\n');
fprintf(fid2,'run 50000');
fprintf(fid2,'\n');
fprintf(fid2,'unfix 1');
fprintf(fid2,'\n');
fprintf(fid2,'\n');
fprintf(fid2,'fix 1 all npt temp $T $T 0.1 x 0.0 0.0 1.0 y 0.0 0.0 1.0 couple xy');
fprintf(fid2,'\n');
fprintf(fid2,'run 50000');
fprintf(fid2,'\n');
fprintf(fid2,'unfix 1');
fprintf(fid2,'\n');
fprintf(fid2,'unfix 2');
fprintf(fid2,'\n\n');
fprintf(fid2,'fix NVE all nve\n');
fprintf(fid2,'fix				ts all temp/rescale 1 $T $T $(100.0*dt) 1.0\n');
fprintf(fid2,'run 50000\n');
fprintf(fid2,'unfix     ts\n');
fprintf(fid2,'unfix     NVE\n\n');

fprintf(fid2,'# ======= write output ======= #\n\n');

fprintf(fid2,'write_data post_equil.data');
fprintf(fid2,'\n');
fclose(fid2);
end

