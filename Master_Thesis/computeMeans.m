function [cstR0,cstY0,cstI0,cstP0,csts0,cstf0,cstW0,cstdf0,cstds0,csts600] = computeMeans(sd,fd)

load estimationData;

smpl   = sd:fd;
cstR0   = mean(R_o(smpl));
cstY0   = mean(dY_o(smpl));
cstI0   = mean(dI_o(smpl));
cstP0   = mean(dP_o(smpl));
csts0   = mean(s_o(smpl));
cstf0   = mean(f_o(smpl));
cstW0   = mean(dW_o(smpl));
cstds0  = mean(ds_o(smpl));
cstdf0  = mean(df_o(smpl));
csts600 = mean(sh60_o(smpl));
