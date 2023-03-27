% Scripts for looking at BCGS_IRO_LS
addpath(genpath('../'));
rng(4)

XXdim = [1000 120 2];
skel = {'bcgs_pip', 'bmgs', 'bcgs_iro', 'bcgs_iro_ls'};
musc = {'cholqr', 'cgs', 'mgs', 'houseqr'};

disp('First observe block behavior')
LaeuchliBlockKappaPlot(XXdim, [], skel, musc);
MonomialBlockKappaPlot(XXdim, [], skel, musc);

%% Set up additional, difficult test matrices