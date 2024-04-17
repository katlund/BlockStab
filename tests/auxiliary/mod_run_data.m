function run_data = mod_run_data(run_data, ind)
% run_data = MOD_RUN_DATA(run_data, ind) is a subroutine for extracting a
% subset of run_data entries according to integer indices ind for
% replotting with GEN_PLOT.
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
run_data.loss_ortho = run_data.loss_ortho(:,ind);
run_data.rel_res = run_data.rel_res(:,ind);
run_data.rel_chol_res = run_data.rel_chol_res(:,ind);
run_data.lgd = run_data.lgd(ind);
run_data.alg_cmap = run_data.alg_cmap(ind,:);
run_data.symb = run_data.symb(ind);
end