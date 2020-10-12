
function prepare_data_Cyprus(Cyprus_data,data_all)
% this function prepares Cyprus data and adds it to the main data file

% read all the sheets
filename =  Cyprus_data;

[~,sheets] = xlsfinfo(filename);
N = length(sheets) ;  
file = cell(N,1) ;
for i = 1:N
    [num,year,ro] = xlsread(filename,i);
    file{i} = num ;    
end
[~,year,~] = xlsread(filename,1);

% output, employment, household expenditures, yhat
n = file{4,1};
y = file{1,1};
l = file{2,1};
x = file{3,1};
yhat = (y(3,1) * n(3,2) * y(:,2) .* n(:,1)) ./ (n(3,1) * y(3,2) * n(:,2));
data = [];
data(:,1) = y(:,1) ./ yhat;
data(:,2) = l(:,1) ./ ((l(3,1) * n(3,2) * l(:,2) .* n(:,1)) ./ (n(3,1) * l(3,2) * n(:,2)));
data(:,3) = x(:,1) ./ ((x(3,1) * n(3,2) * x(:,2) .* n(:,1)) ./ (n(3,1) * x(3,2) * n(:,2)));
data(:,4) = yhat;

% household debt, government debt, government expenditure, government
% transfers, net export, current account, net foreign assets, government
% interest payments, household interest payments, government tax revenue,
% recapitalization costs, gross export
for i = 5:16
   data(:,i) = file{i,1}(:,1) ./ yhat;
end

% correct government debt
for i = 2:13
    data(i,6) = data(i-1,6) + data(i,7) + data(i,8) - data(i,14);
end


% RHO

% long_term_rates speads
for i =1:length(file{17,1})
    ltr_median(i) = median(file{17,1}(:,i));
    file{17,1}(:,i) = file{17,1}(:,i) - ltr_median(i);
end
Cyprus = file{17,1}(12,:)/100;
ltr = file{17,1}/100;
ltr1_add1 = ltr(1:11,1:8);
ltr1_add1 = ltr1_add1(:);
ltr1_add2 = ltr(:,9:13);
ltr1_add2 = ltr1_add2(:);
ltr1(:,2) = [ltr1_add1; ltr1_add2];
ltr1(:,1) = ones(length(ltr1),1);
ltr2_add1 = ltr(1:11,1:5);
ltr2_add1 = ltr2_add1(:);
ltr2_add2 = ltr(:,6:10);
ltr2_add2 = ltr2_add2(:);
ltr2(:,2) = [ltr2_add1; ltr2_add2];
ltr2(:,1) = ones(length(ltr2),1);

% SME loans rates dpreads 
for i = 1:8
    ltr_median(i) = median(file{18,1}(1:11,i));
    sme1(:,i) = file{18,1}(1:11,i) - ltr_median(i);
end
for i = 9:13
    ltr_median(i) = median(file{18,1}(:,i));
    sme2(:,i) = file{18,1}(:,i) - ltr_median(i);
end
sme2 = sme2(:,9:13);
sme = [sme1(:); sme2(:)]/100;

% deposit rates spreads
for i = 1:5
    ltr_median(i) = median(file{19,1}(1:11,i));
    dr1(:,i) = file{19,1}(1:11,i) - ltr_median(i);
end
for i = 6:10
    ltr_median(i) = median(file{19,1}(:,i));
    dr2(:,i) = file{19,1}(:,i) - ltr_median(i);
end
dr2 = dr2(:,6:10);
dr = [dr1(:); dr2(:)]/100;

% regression
reg{1,1} = regress(sme,ltr1);
reg{2,1} = regress(dr,ltr2);
for j = 1:2
    for i = 1:length(Cyprus)
        rates(i,j) = reg{j,1}(1,1) + reg{j,1}(2,1)*Cyprus(i);
    end
end
data(:,17) = mean(rates,2);


% identification of bhat with the US states
data(:,18) = data(:,5);
debt = file{23,1} ./ file{24,1} ./ 1000000;
coef9 = regress(debt(:,7), [debt(:,6), debt(:,5), debt(:,3), debt(:,1)]);
data(10,18) = coef9(1)*data(9,6) + coef9(2)*data(8,6) + coef9(3)*data(6,6) + coef9(4)*data(4,6);
coef10 = regress(debt(:,8), [debt(:,7), debt(:,5), debt(:,3), debt(:,1)]);
data(11,18) = coef10(1)*data(10,6) + coef10(2)*data(8,6) + coef10(3)*data(6,6) + coef10(4)*data(4,6);
coef11 = regress(debt(:,9), [debt(:,8), debt(:,5), debt(:,3), debt(:,1)]);
data(12,18) = coef11(1)*data(11,6) + coef11(2)*data(8,6) + coef11(3)*data(6,6) + coef11(4)*data(4,6);
coef12 = regress(debt(:,10), [debt(:,9), debt(:,5), debt(:,3), debt(:,1)]);
data(13,18) = coef12(1)*data(12,6) + coef12(2)*data(8,6) + coef12(3)*data(6,6) + coef12(4)*data(4,6);

% wage and price
for i = 20:22
    file{i,1}(13,:) = sum(file{i,1}(1:12,:))/12;
    file{i,1}(1:12,:) = file{i,1}(1:12,:) ./ file{i,1}(13,:);
    file{i,1} = file{i,1} ./ file{i,1}(:,3);
    file{i,1} = file{i,1}(1:12,:);
end
data(:,19) = ((file{22,1}(8,:) + file{20,1}(8,:) ./ file{21,1}(8,:))/2)';
data(:,20) = 0;

% merge with the main dataset
data(:,2:end+1) = data(:,1:end);
data(:,1) = string(year(2:14,1)');
base = xlsread(data_all);
fin = [base;data];
fin(144:156,21) = fin(131:143,21);
writematrix(fin,'data.xlsx')
end

