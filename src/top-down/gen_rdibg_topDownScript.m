%% Afnan Mostafa
%% 03/16/2022
%% code to generate randomly-distributed interlayer-bonded bilayer graphene (rd-ibg)
%% from diamane (2D diamond) structures %%%

%% %%% clear variables %%% %%
clear
clc
close all

%% %%% reading the data file and storing spatial coordinate values %%% %%
layers = 2;
stackings = "aa";
lengths = 10;
widths = 10;
stackings = lower(stackings);

%% reading the data file and storing spatial coordinate values
for h=1:length(lengths)
    for g=1:length(widths)
        for idx=1:length(stackings)
            %% edge case 1:
            if sum(strcmpi(stackings(idx), "ab") || strcmpi(stackings(idx), "aa")) ~= 1
                error('wrong stacking, use either aa or ab')
            end
            %% edge case 2:
            if layers ~= 2
                error('use 2 layers')
            end
            %% read atom data
            length_sheet = lengths(h);
            width_sheet = widths(g);
            [file,len,wid,~] = generate_graphene(length_sheet,width_sheet,stackings(idx),layers);
            [s] = readtextfile(file,5,15,'','#');
            datafile = file;
            filename = file;
            index=s{1};
            atom_type=s{2};
            x = s{3};
            y = s{4};
            z = s{5};
            zz_atom_oneStrip = ceil((wid*10)/2.45);
            ac_atom_oneStrip = ceil((len*10)/4.2);
            totAtom_zz = zz_atom_oneStrip*4;
            totAtom_ac = ac_atom_oneStrip*4;
            top_end = length(x);
            top_strt = length(x) - totAtom_zz;
            bot_strt = 1;
            bot_end = bot_strt + totAtom_zz;
            
            %% %%% selecting the aligned atoms %%% %%
            x_bottom = x(1:length(x)/2);
            x_top = x(length(x)/2+1:end);
            
            y_bottom = y(1:length(y)/2);
            y_top = y(length(y)/2+1:end);
            
            z_bottom = z(1:length(z)/2);
            z_top = z(length(z)/2+1:end);
            
            bottom_aligned_atoms = zeros(length(x_top),1);
            top_aligned_atoms = zeros(length(x_top),1);
            
            counter1 = 0;
            for i = 1:length(x_bottom)
                for j = 1:length(x_top)
                    if (y_bottom(i) == y_top(j)) && (x_bottom(i) == x_top(j))
                        counter1 = counter1 + 1;
                        bottom_aligned_atoms(counter1,1) = index(i);
                        top_aligned_atoms(counter1,1) = index(j);
                    end
                end
            end
            
            bottom_aligned_atoms = bottom_aligned_atoms(1:counter1);
            top_aligned_atoms = top_aligned_atoms(1:counter1);
            top_aligned_atoms_added = top_aligned_atoms + length(x_bottom);
            
            id_cc_b = index((bottom_aligned_atoms));
            id_cc_t = index(top_aligned_atoms_added);
            
            %% %%% indices of atoms for hydrogenation %%% %%
            if strcmpi(stackings(idx), "ab")
                %% ################ AB ############ %%
                %% bottom layer %%
                A = index(1:length(x_bottom));
                set = (id_cc_b);
                %% top layer %%
                A1 = index(length(x_top)+1:end);
                set1 = (id_cc_t);
                bound = index((y<min(y)+1 | x<min(x)+0.5 | x>max(x)-0.5 | y>max(x)-0.5));
            elseif strcmpi(stackings(idx), "aa")
                %% ################ AA ############ %%
                for u = 1:length(atom_type)
                    if (rem((index(u,1)), 2) == 0)
                        A(u,1) = u;
                    else
                        A1(u,1) = index(u,1);
                    end
                end
                
                A(A==0) = [];
                A1(A1==0) = [];
                ch_bot=A1(A1 < (length(x_bottom)+1));
                ch_top=A1(A1 > length(x_top));
                set=A(A < (length(x_bottom)+1));
                set1=A(A > (length(x_bottom)));
                bound = index((y<min(y)+1 | x<min(x)+0.5 | x>max(x)-0.5 | y>max(x)-0.5));
                
            end
            
            %% top-down from diamane to rd-ibg
            bond_density = 1.25;
            bonddens_remove = 50-bond_density;
            bonddens_remove = bonddens_remove*2;
            bonddens_frac = bonddens_remove/100;
            final_bd = bond_density;
            if strcmpi(stackings(idx), "ab")
                tot_rand_bonds = round(length(set)*bonddens_frac);
                id_rand_bonds = randperm(numel(set),tot_rand_bonds)';
                removed_bot = set(id_rand_bonds);
                removed_top = set1(id_rand_bonds);
                set = set(~ismember(set,removed_bot));
                set1 = set1(~ismember(set1,removed_top));
            elseif strcmpi(stackings(idx), "aa")
                tot_rand_bonds = round(length(set)*bonddens_frac);
                id_rand_bonds = randperm(numel(set),tot_rand_bonds)';
                removed_bot = set(id_rand_bonds(1:end));
                removed_top = set1(id_rand_bonds(1:end));
                set = set(~ismember(set,removed_bot));
                set1 = set1(~ismember(set1,removed_top));
            end
            
            %% %%% Section 2.2: select atoms within range %%% %%
            set = set(~ismember(set,bound));
            set1 = set1(~ismember(set1,bound));
            mid_atom_list = sort(set);
            hydrogented_C_list=zeros(1,1);
            mid_x = x_bottom(mid_atom_list);
            mid_y = y_bottom(mid_atom_list);
            rad = 1.5;
            
            cout = 0;
            for nm = 1: length(mid_atom_list)
                for ia=1:length(x_bottom)
                    if mid_atom_list(nm) ~= ia
                        if (((x_bottom(ia)-mid_x(nm))^2)+(y_bottom(ia)-mid_y(nm))^2) < rad^2
                            cout = cout+1;
                            hydrogented_C_list(cout,1) = ia;
                        end
                    end
                end
            end
            hydrogented_C_list = unique(hydrogented_C_list);
            
            mid_atom_list_top = set1;
            hydrogented_C_list_top=zeros(1,1);
            mid_x_top = x(mid_atom_list_top);
            mid_y_top = y(mid_atom_list_top);
            cout2 = 0;
            for vv = 1: length(mid_atom_list_top)
                for ib=length(x_top)+1:length(x)
                    if mid_atom_list_top(vv) ~= ib
                        if (((x(ib)-mid_x_top(vv))^2)+(y(ib)-mid_y_top(vv))^2) < rad^2
                            cout2 = cout2+1;
                            hydrogented_C_list_top(cout2,1) = ib;
                        end
                    end
                end
            end
            hydrogented_C_list_top = unique(hydrogented_C_list_top);
            
            %% %%% writing LAMMPS input script for C-C atoms to be moved towards each other in order to be bonded %%% %%
            finalfilename = 'RDIBG_%s_%dlayers_%dnmx%dnm_%.2fBD.data';
            finalfile = sprintf(finalfilename,stackings(idx),layers,length_sheet,width_sheet,final_bd);
            
            script = 'in_%s_%dnmx%dnm_%0.2fBD.lmp';
            fid2 = sprintf(script, stackings(idx),length_sheet,width_sheet,bond_density);
            fid2 = fopen(fid2, 'w');
            fprintf(fid2, '########### Prepare: RD-IBG (%0.2f%% bond density, %s) from 2D diamond (top-down) ###########',bond_density,stackings);
            fprintf(fid2,'\n\n');
            
            fprintf(fid2,'###============= simulation setup =============###\n\n');
            fprintf(fid2,'units         metal\n');
            fprintf(fid2,'dimension     3\n');
            fprintf(fid2,'boundary      p p p\n');
            fprintf(fid2,'processors    * * 1\n\n');
            
            fprintf(fid2,'###============= neighbor list build-up settings =============###\n\n');
            fprintf(fid2,'neighbor      2.0 bin\n');
            fprintf(fid2,'neigh_modify  every 1 delay 0 check yes\n\n');
            
            fprintf(fid2,'###============= timestep size in ps =============###\n\n');
            fprintf(fid2,'timestep      0.0001\n');
            fprintf(fid2,'variable      dt equal 0.0001\n\n');
            
            fprintf(fid2,'###============= temp and seed variables =============###\n\n');
            fprintf(fid2,'variable      ini_temp equal 10.0\n');
            fprintf(fid2,'variable      eq_temp equal 300.0\n');
            fprintf(fid2,'variable      st_temp equal 300.0\n');
            fprintf(fid2,'variable      velSeed equal 42618127\n\n');
            
            fprintf(fid2,'###============= data file =============###\n\n');
            fprintf(fid2,'atom_style	atomic\n');
            fprintf(fid2,'read_data       %s',finalfile);
            fprintf(fid2,'\n\n');
            fprintf(fid2,'###============= potential file information =============###\n\n');
            fprintf(fid2,'pair_style	airebo 3 1 1\n');
            fprintf(fid2,'pair_coeff 	* * CH.airebo C C H\n\n');
            
            fprintf(fid2,'###============= stress characterization =============###\n\n');
            fprintf(fid2,'group         hydro type 3\n');
            fprintf(fid2,'group         carbon type 1 2\n\n');
            fprintf(fid2,'compute       1 all stress/atom NULL\n');
            fprintf(fid2,'compute       2 all reduce sum c_1[1] c_1[2] c_1[4]\n\n');
            fprintf(fid2,'compute       3 all displace/atom\n');
            fprintf(fid2,'compute       4 carbon coord/atom cutoff 1.7 group all\n');
            fprintf(fid2,'compute       5 hydro coord/atom cutoff 1.3 group carbon \n\n');
            fprintf(fid2,'variable      Devery equal 10\n');
            fprintf(fid2,'variable      Drepeat equal 500\n');
            fprintf(fid2,'variable      Dfreq equal "(v_Devery*v_Drepeat)"\n\n');
            
            fprintf(fid2,'###============= dump =============###\n\n');
            fprintf(fid2,'variable      svm atom sqrt((((c_1[1]-c_1[2])^2+(c_1[2]-c_1[3])^2+(c_1[3]-c_1[1])^2)+6*(c_1[4]^2+c_1[5]^2+c_1[6]^2))/2)\n\n');
            fprintf(fid2,'fix           avgs1 all ave/atom ${Devery} ${Drepeat} ${Dfreq} v_svm c_1[*] c_3[*] c_4 c_5\n\n');
            fprintf(fid2,'dump          mini all custom ${Dfreq} minimize.xyz id type x y z v_svm f_avgs1[*]\n');
            fprintf(fid2,'dump_modify   mini sort id\n\n');
            
            fprintf(fid2,'group         atom1 id \t');
            
            if strcmpi(stackings(idx), "ab")
                fprintf(fid2, '%d ',set);
                fprintf(fid2,'\n\n');
                fprintf(fid2,'group         atom2 id \t');
                fprintf(fid2, '%d ',set1);
                
            elseif strcmpi(stackings(idx), "aa")
                fprintf(fid2, '%d ',set);
                fprintf(fid2,'\n\n');
                fprintf(fid2,'group         atom2 id \t');
                fprintf(fid2, '%d ',set1);
            end
            
            fprintf(fid2,'\n\n');
            fprintf(fid2,'thermo        100\n\n');
            fprintf(fid2,'displace_atoms atom1 move 0 0 1.15 units box\n');
            fprintf(fid2,'displace_atoms atom2 move 0 0 -1.15 units box\n\n');
            fprintf(fid2,'run 1000\n\n');
            
            fprintf(fid2,'thermo_style  custom step time temp press pxx pyy pzz pe ke etotal\n');
            fprintf(fid2,'thermo        1000\n\n');
            
            fprintf(fid2,'###============= minimization =============###\n\n');
            fprintf(fid2,'min_style		cg\n');
            fprintf(fid2,'minimize		0 1e-5 200000 200000\n\n');
            
            fprintf(fid2,'min_style		cg\n');
            fprintf(fid2,'fix 			1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz\n');
            fprintf(fid2,'minimize		0 1e-8 200000 200000\n');
            fprintf(fid2,'unfix			1\n\n');
            
            fprintf(fid2,'write_data    post_cg1.data\n\n');
            
            fprintf(fid2,'write_dump    all atom dump1.atom\n');
            fprintf(fid2,'displace_atoms all random 0.05 0.05 0.05 91897\n');
            fprintf(fid2,'write_dump    all atom dump2.atom\n\n');
            
            fprintf(fid2,'min_style		cg\n');
            fprintf(fid2,'fix 			1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz\n');
            fprintf(fid2,'minimize		0 1e-8 200000 200000\n');
            fprintf(fid2,'unfix			1\n\n');
            
            fprintf(fid2,'write_data    post_cg2.data\n\n');
            fprintf(fid2,'write_restart minimized_struc.restart\n\n');
            
            fprintf(fid2,'###============= velocity distribution =============###\n\n');
            fprintf(fid2,'velocity      all create ${ini_temp} ${velSeed} mom yes rot yes dist gaussian units box\n\n');
            fprintf(fid2,'###============= equilibration =============###\n\n');
            fprintf(fid2,'thermo        1000\n\n');
            
            fprintf(fid2,'fix           NPT all npt temp ${ini_temp} ${ini_temp} 0.01 x 0.0 0.0 0.1 y 0.0 0.0 0.1 drag 2.0\n');
            fprintf(fid2,'run           100000\n');
            fprintf(fid2,'unfix         NPT\n\n');
            
            fprintf(fid2,'fix           NPT all npt temp ${ini_temp} ${eq_temp} 0.01 x 0.0 0.0 0.5 y 0.0 0.0 0.5 drag 2.0\n');
            fprintf(fid2,'run           120000\n');
            fprintf(fid2,'unfix         NPT\n\n');
            
            fprintf(fid2,'fix           NPT all npt temp ${eq_temp} ${eq_temp} 0.01 x 0.0 0.0 0.1 y 0.0 0.0 0.1 drag 2.0\n');
            fprintf(fid2,'run           100000\n');
            fprintf(fid2,'unfix         NPT\n\n');
            
            fprintf(fid2,'fix           1 all nvt temp ${st_temp} ${st_temp} $(100.0*dt) drag 2.0\n');
            fprintf(fid2,'run           100000\n');
            fprintf(fid2,'unfix         1\n\n');
            
            fprintf(fid2,'undump        mini\n');
            fprintf(fid2,'unfix         avgs1\n\n');
            
            fprintf(fid2,'write_data    post_equil.data\n\n');
            fprintf(fid2,'write_restart relaxed_struc.restart\n\n');
            fprintf(fid2,'print "Minimization + Equilibration DONE"\n');
            fclose(fid2);
            
            %% %%% data file including the H atoms (No C atoms in this data file) %%% %%
            output = 'RDIBG_%s_%dx%d_%.2fBD_onlyH.data';
            outfile = sprintf(output,stackings(idx),length_sheet,width_sheet,final_bd);
            fid = fopen(outfile,'w');
            
            if strcmpi(stackings(idx), "ab")
                l=1;
                cout = length(hydrogented_C_list);
                for m=1:cout
                    fprintf(fid,'%g 3 %g %g %g\n',length(atom_type)+l,x_bottom(hydrogented_C_list(m)), y_bottom(hydrogented_C_list(m)),z_bottom(m)-0.9);
                    l=l+1;
                end
                clear m l
                l = 1;
                cout2 = length(hydrogented_C_list_top);
                for m=1:cout2
                    fprintf(fid,'%g 3 %g %g %g\n',length(atom_type)+cout+l,x(hydrogented_C_list_top(m)), y(hydrogented_C_list_top(m)),z_top(m)+0.9);
                    l=l+1;
                end
                fclose(fid);
                total_atoms = length(index)+length(hydrogented_C_list)+length(hydrogented_C_list_top);
                types = 3;
            elseif strcmpi(stackings(idx), "aa")
                l=1;
                cout = length(hydrogented_C_list);
                for m=1:cout
                    fprintf(fid,'%g 3 %g %g %g\n',length(atom_type)+l,x_bottom(hydrogented_C_list(m)), y_bottom(hydrogented_C_list(m)),z_bottom(hydrogented_C_list(m))-0.9);
                    l=l+1;
                end
                clear l m
                l=length(hydrogented_C_list_top)+1;
                cout2 = length(hydrogented_C_list_top);
                for m=1:cout2
                    fprintf(fid,'%g 3 %g %g %g\n',length(atom_type)+l,x(hydrogented_C_list_top(m)), y(hydrogented_C_list_top(m)),z(hydrogented_C_list_top(m))+0.9);
                    l=l+1;
                end
                fclose(fid);
                total_atoms = length(index)+length(hydrogented_C_list)+length(hydrogented_C_list_top);
                types = 3;
            end
            
            %% final modification to output file
            onetimerun= true;
            if onetimerun == true
                str5 = sprintf("rm %s '",finalfile);
                command5 = (str5);
                system(command5);
                
                command = 'cat %s %s > %s';
                command = sprintf(command,filename,outfile, finalfile);
                system(command);
                
                command = 'sed -i "s/2 atom types/3 atom types/" %s';
                command = sprintf(command, finalfile);
                system(command);
                
                command = 'sed -i "s/%d atoms/%d atoms/1" %s';
                command = sprintf(command, length(index), total_atoms, finalfile);
                system(command);
                
                command = 'sed -i "13 a 3 1.00784" %s';
                command = sprintf(command, length(index), total_atoms, finalfile);
                system(command);
                
                % comment the following 3 lines if you want to see the H
                % data in a separate file
                str = sprintf("rm %s %s'",outfile,file);
                command = (str);
                system(command);
                
                clearvars -except g h idx stackings lengths widths layers
            end
        end
    end
end

fclose('all');
%% %%%                                            End of script                                                             %%% %%