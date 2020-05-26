
%% Review Questions

% Trade balance is the difference between country's imports and exports for
% some period.

% Current account is country's trade balance plus net income and direct payments.

% Stylized business cycle facts:
% 1. High global volatility: the cross-country average standard deviation
% of output is twice as large as its U.S. counterpart.
% 2. Excess consumption volatility: private consumption (including
% durables) tends to be more volatile than output.
% 3. Global ranking of volatilities: the ranking of cross-country average standard deviations from top to bottom is imports, investment, exports, government spending, consumption, and output.
% 4. Cyclicality: consumption, investment, exports, and imports are
% procyclical.
% 5. Cyclicality: the trade balance and the current account are
% countercyclical.
% 6. Cyclicality: the share of government consumption in output is acyclical.
% 7. Persistence: All components of demand (c, g, i, x) and supply (y, m) are positively serially correlated.
% 8. Excess volatility if poor or emerging: Business Cycles in poor and
% emerging countries are about twice as volatile as in rich countries.
% 9. Government spending more countercyclical if richer: the share of government consumption is countercyclical in rich countries, but acyclical in emerging and poor countries.


%% Exercises. Country - Russia, period: 1994-2017
clear all
clc
close all

%% read the data 
[data, columns] = xlsread('data.xlsx');

var_init =  NaN([length(data) 6]);

% log of deflated output per capita
var_init(:,1) = log(data(:,4) / data(:,10) * data(:,5));
% log of deflated consumption per capita
var_init(:,2) = log(data(:,8) / data(:,10) * data(:,5));
% log of deflated government consumption per capita
var_init(:,3) = log(data(:,6) / data(:,10) * data(:,5));
% log of deflated investment per capita
var_init(:,4) = log(data(:,7) / data(:,10) * data(:,5));
% log of deflated export per capita
var_init(:,5) = log(data(:,3) / data(:,10) * data(:,5));
% log of deflated import per capita
var_init(:,6) = log(data(:,9) / data(:,10) * data(:,5));
% time array
t = [data(:,1), data(:,1).^2];

%% detrend and find percetage deviations

var_cyc = NaN([length(data) 6]);
var_sec = NaN([length(data) 6]);
var_fin = NaN([length(data) 6]);

for i = 1:6
    [a,b,var_cyc(:,i)] = regress(var_init(:,i),t);
    var_sec(:,i) = var_init(:,i)-var_cyc(:,i);
    var_fin(:,i) = var_cyc(:,i)./var_sec(:,i);
end

%% find trade balance and current account
var_init(:,7) = (data(:,3) / data(:,10) * data(:,5) - data(:,9) / data(:,10) * data(:,5))./exp(var_sec(:,1));
[a,b,var_fin(:,7)] = regress(var_init(:,7),t);

var_init(:,8) = ((data(:,2)/100.*data(:,4))/ data(:,10) * data(:,5))./exp(var_sec(:,1));
[a,b,var_fin(:,8)] = regress(var_init(:,8),t);

%% display necessary results
disp('STYLIZED FACTS FOR RUSSIA, PERIOD 1994-2017');
fprintf(1, '\n');
disp('Fact 1. Average standard deviation of output is a bit larger than in the US:');
fprintf('sigma_y: %0.1f %%', std(var_fin(:,1))*100);
fprintf(1, '\n');
fprintf(1, '\n');
disp('Fact 2. Private consumption is a bit more violate than the output:')
fprintf('sigma_c/sigma_y: %0.3f %', std(var_fin(:,2))/std(var_fin(:,1)));
fprintf(1, '\n');
fprintf(1, '\n');
disp('Fact 3. The ranking of cross-country average standard deviations from top to bottom is exports, imports, investment, government spending, consumption, and output:')
fprintf('sigma_x/sigma_y: %0.1f %', std(var_fin(:,5))/std(var_fin(:,1)));
fprintf(1, '\n');
fprintf('sigma_m/sigma_y: %0.1f %', std(var_fin(:,6))/std(var_fin(:,1)));
fprintf(1, '\n');
fprintf('sigma_i/sigma_y: %0.1f %', std(var_fin(:,4))/std(var_fin(:,1)));
fprintf(1, '\n');
fprintf('sigma_g/sigma_y: %0.1f %', std(var_fin(:,3))/std(var_fin(:,1)));
fprintf(1, '\n');
fprintf('sigma_c/sigma_y: %0.1f %', std(var_fin(:,2))/std(var_fin(:,1)));
fprintf(1, '\n');
fprintf(1, '\n');
disp('Facts 4-6. All variables are procyclical:')
fprintf('corr(c,y): %0.2f %', corr(var_fin(:,2),var_fin(:,1)));
fprintf(1, '\n');
fprintf('corr(i,y): %0.2f %', corr(var_fin(:,4),var_fin(:,1)));
fprintf(1, '\n');
fprintf('corr(x,y): %0.2f %', corr(var_fin(:,5),var_fin(:,1)));
fprintf(1, '\n');
fprintf('corr(m,y): %0.2f %', corr(var_fin(:,6),var_fin(:,1)));
fprintf(1, '\n');
fprintf('corr(tb,y): %0.2f %', corr(var_fin(:,7),var_fin(:,1)));
fprintf(1, '\n');
fprintf('corr(ca,y): %0.2f %', corr(var_fin(:,8),var_fin(:,1)));
fprintf(1, '\n');
fprintf('corr(g,y): %0.2f %', corr(var_fin(:,3),var_fin(:,1)));
fprintf(1, '\n');
fprintf(1, '\n');
disp('Fact 7. All components of demand (c, g, i, x) and supply (y, m) are positively serially correlated:')
fprintf('corr(yt,yt-1): %0.2f %', corr(var_fin(2:24,1),var_fin(1:23,1)));
fprintf(1, '\n');
fprintf('corr(ct,ct-1): %0.2f %', corr(var_fin(2:24,2),var_fin(1:23,2)));
fprintf(1, '\n');
fprintf('corr(gt,gt-1): %0.2f %', corr(var_fin(2:24,3),var_fin(1:23,3)));
fprintf(1, '\n');
fprintf('corr(it,it-1): %0.2f %', corr(var_fin(2:24,4),var_fin(1:23,4)));
fprintf(1, '\n');
fprintf('corr(xt,xt-1): %0.2f %', corr(var_fin(2:24,5),var_fin(1:23,5)));
fprintf(1, '\n');
fprintf('corr(mt,mt-1): %0.2f %', corr(var_fin(2:24,6),var_fin(1:23,6)));
fprintf(1, '\n');