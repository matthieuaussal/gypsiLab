function Mh = hmxSingle(Mh)
%+========================================================================+
%|                                                                        |
%|               OPENHMX, H-MATRIX COMPRESSION AND ALGEBRA                |
%|              openHmx is part of GYPSYLAB toolbox - v0.20               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
%| Licence Creative Commons BY-NC-SA 4.0, Attribution, NonCommercial and  |
%| ShareAlike (see http://creativecommons.org/licenses/by-nc-sa/4.0/).    |
%| This software is the property from Centre de Mathematiques Appliquees  |
%| de l'Ecole polytechnique, route de Saclay, 91128 Palaiseau, France.    |
%|                                                            _   _   _   |
%| Please acknowledge the GYPSILAB toolbox in programs       | | | | | |  |
%| or publications in which you use the code. For openHmx,    \ \| |/ /   |
%| we suggest as reference :                                   \ | | /    |
%| [1] : www.cmap.polytechnique.fr/~aussal/gypsilab             \   /     |
%| [2] : 13th International Conference on Mathematical           | |      |
%| and Numerical Aspects of Wave Propagation, University of      | |      |
%| Minnesota, may 2017. "OpenHmX, an open-source H-Matrix        | |      |
%| toolbox in Matlab".                                           | |      |
%|_______________________________________________________________|_|______|
%| Author(s)  : Matthieu Aussal - CMAP, Ecole polytechnique               |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Convert H-Matrix to single precision                      |
%+========================================================================+

% H-Matrix (recursion)
if (Mh.typ == 0)
    for i = 1:4
        Mh.chd{i} = hmxSingle(Mh.chd{i});
        Mh.row{i} = single(Mh.row{i});
        Mh.col{i} = single(Mh.col{i});
    end
    
% Compressed leaf
elseif (Mh.typ == 1)
    Mh.dat{1} = single(Mh.dat{1});
    Mh.dat{2} = single(Mh.dat{2});
    
% Full leaf
elseif (Mh.typ == 2)
    Mh.dat = single(Mh.dat);
    
% Sparse leaf
elseif (Mh.typ == 3)
    Mh.dat = Mh.dat;
    
% Unknown type
else
    error('hmxSingle.m : unavailable case')
end

end