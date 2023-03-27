% Scripts for looking at BCGS_IRO_LS
addpath(genpath('../'));
rng(4)

XXdim = [1000 120 2];
skel = {'bcgs_pip', 'bmgs', 'bcgs_iro', 'bcgs_iro_ls', 'bcgs_iro_bl'};
musc = {'houseqr'};

disp('First observe block behavior')
LaeuchliBlockKappaPlot(XXdim, [], skel, musc);
MonomialBlockKappaPlot(XXdim, [], skel, musc);