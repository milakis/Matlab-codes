clear all, close all, clc
addpath('/usr/lib/dynare/matlab') ;

%closed economy
dynare Smoother_closed
irfs_cl=cell2mat(struct2cell(oo_.irfs));
irfs_cl_unins = irfs_cl([2 8 9 12 13 15 17 18 19 20 21 25 27 30],:);
irfs_cl_ins = irfs_cl([37,43,44,47,48,50,52,53,54,55,56,60,62,65],:);
T=20;
which_plot_ = strvcat('Interest rate', 'Firm owners consumption', 'Employed workers consumption', 'Aggregate consumption','Firm owners assets', 'Employed workers total assets', 'Matches', 'Vacancies', 'Job-finding rate',' Job-loss rate','Employment', 'Output', 'Labor servies price', 'Inflation'); 
nplot       = size(which_plot_,1);
nrow        = 4;            
ncol = ceil(nplot/nrow); 
figure('Name','Closed economy');
for jkl=1:nplot
    subplot(nrow,ncol,jkl);
    hold on;
    plot(1:T,irfs_cl_unins(jkl,1:T));
    hold on;
    plot(1:T,irfs_cl_ins(jkl,1:T));
    title(deblank(which_plot_(jkl,:)));
end
legend('Imperfect insurance','Perfect insurance');

%open economy with flexible exchange rate
dynare Smoother_open_flex
irfs_op=cell2mat(struct2cell(oo_.irfs));
irfs_op_unins = irfs_op([2 8 9 12 13 15 17 18 19 20 21 25 27 30 36 37 38 39],:);
irfs_op_ins = irfs_op([41,47,48,51,52,54,56,57,58,59,60,64,66,69,75,76,77,78],:);
T=20;
which_plot_ = strvcat('Interest rate', 'Firm owners consumption', 'Employed workers consumption', 'Aggregate consumption','Firm owners assets', 'Employed workers total assets', 'Matches', 'Vacancies', 'Job-finding rate',' Job-loss rate','Employment', 'Output', 'Labor servies price', 'Inflation', 'Trade', 'Exchange rate','NFA','Net exports'); 
nplot       = size(which_plot_,1);
nrow        = 5;            
ncol = ceil(nplot/nrow); 
figure('Name','Open economy with flexible exchange rate');
for jkl=1:nplot
    subplot(nrow,ncol,jkl);
    hold on;
    plot(1:T,irfs_op_unins(jkl,1:T));
    hold on;
    plot(1:T,irfs_op_ins(jkl,1:T));
    title(deblank(which_plot_(jkl,:)));
end
legend('Imperfect insurance','Perfect insurance');

%open economy with fixed exchange rate
dynare Smoother_open_peg
irfs_op_peg=cell2mat(struct2cell(oo_.irfs));
copy=repmat(0,1,100);
irfs_op_peg = [irfs_op_peg(1,:) ; copy ; irfs_op_peg(2:35,:) ; copy ; irfs_op_peg(36:38,:) ; copy ; irfs_op_peg(39:72,:) ; copy ; irfs_op_peg(73:end,:)];
irfs_op_peg_unins = irfs_op_peg([2 8 9 12 13 15 17 18 19 20 21 25 27 30 36 37 38 39],:);
irfs_op_peg_ins = irfs_op_peg([41,47,48,51,52,54,56,57,58,59,60,64,66,69,75,76,77,78],:);
figure('Name','Open economy with currency peg');
for jkl=1:nplot
    subplot(nrow,ncol,jkl);
    hold on;
    plot(1:T,irfs_op_peg_unins(jkl,1:T));
    hold on;
    plot(1:T,irfs_op_peg_ins(jkl,1:T));
    title(deblank(which_plot_(jkl,:)));
end
legend('Imperfect insurance','Perfect insurance');

%calculate difference
dif_closed = (irfs_cl_ins - irfs_cl_unins)./irfs_cl_ins;
irfs_op_peg_unins = irfs_op_peg([2 8 9 12 13 15 17 18 19 20 21 25 27 30],:);
irfs_op_peg_ins = irfs_op_peg([41,47,48,51,52,54,56,57,58,59,60,64,66,69],:);
dif_peg = (irfs_op_peg_ins - irfs_op_peg_unins)./irfs_op_peg_ins;

