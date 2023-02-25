function [indx_rand] = gen_rand_bot(bot_align_atoms,bond_density)
%   creates indices for random bonds

bond_den_frac = bond_density/100;
rand_bonds = ceil(length(bot_align_atoms)*bond_den_frac);
indx_rand = randperm(numel(bot_align_atoms),rand_bonds)';

end

