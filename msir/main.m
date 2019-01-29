% main.m
% "/home/rasoul/software/matlab2006a/bin/matlab" -nodesktop

global climate_norm
global xdata
global itter


itter = 0;
cd '/home/rasoul/Dropbox/Programming/iran_paper/code/msir';
ADDRESS = '/home/rasoul/Dropbox/Programming/iran_paper/data.xlsx';
data1 = xlsread(ADDRESS);
INTERVAL = 50:70;
data = data1(INTERVAL,2);
%temp_new = data;

%%%%%%%%%%%%%%%%%%%% Climate Data %%%%%%%%%%%%%%%%%%%%%%
ADRS_SHUM_13 = '/home/rasoul/Dropbox/Programming/iran_paper/data/shum_2013.csv';
ADRS_SHUM_14 = '/home/rasoul/Dropbox/Programming/iran_paper/data/shum_2014.csv';
ADRS_SHUM_15 = '/home/rasoul/Dropbox/Programming/iran_paper/data/shum_2015.csv';

shum_13 = csvread(ADRS_SHUM_13);
shum_14 = csvread(ADRS_SHUM_14);
shum_15 = csvread(ADRS_SHUM_15);

shum = [shum_13(51:53,2); shum_14(1:18,2)];

shum =shum * 1e-5 + .3277;

shum_ave = mean(shum,2);
norm_shum = shum_ave;
m = 2 ./ (max(norm_shum) - min(norm_shum));
norm_shum = m * (norm_shum - min(norm_shum)) - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ADRS_TMP_13 = '/home/rasoul/Dropbox/Programming/iran_paper/data/tmp_2013.csv';
ADRS_TMP_14 = '/home/rasoul/Dropbox/Programming/iran_paper/data/tmp_2014.csv';
ADRS_TMP_15 = '/home/rasoul/Dropbox/Programming/iran_paper/data/tmp_2015.csv';


tmp_13 = csvread(ADRS_TMP_13);
tmp_14 = csvread(ADRS_TMP_14);
tmp_15 = csvread(ADRS_TMP_15);


tmp = [tmp_13(51:53,2); tmp_14(1:18,2)];
tmp = tmp - 273;


tmp_ave = mean(tmp, 2);
norm_tmp = tmp_ave;
m = 2 ./ (max(norm_tmp)- min(norm_tmp));
norm_tmp = m * (norm_tmp - min(norm_tmp)) - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ADRS_PRECIP_13 = '/home/rasoul/Dropbox/Programming/iran_paper/data/precip_2013.csv';
ADRS_PRECIP_14 = '/home/rasoul/Dropbox/Programming/iran_paper/data/precip_2014.csv';
ADRS_PRECIP_15 = '/home/rasoul/Dropbox/Programming/iran_paper/data/precip_2015.csv';

precip_13 = csvread(ADRS_PRECIP_13);
precip_14 = csvread(ADRS_PRECIP_14);
precip_15 = csvread(ADRS_PRECIP_15);

precip = [precip_13(51:53,2); precip_14(1:18,2)];

%precip = precip;

precip_ave = mean(precip, 2);
norm_precip = precip_ave;
m = 2 ./ ( max(norm_precip) - min(norm_precip) );
norm_precip = m * (norm_precip - min(norm_precip)) - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clim_data = [norm_shum norm_tmp norm_precip];

xdata = 1:size(data,1);
ydata = data;
Yfun = @(x, xdata) infection(x, xdata, clim_data);
fun = @(x)sum((Y(x)-ydata).*(Y(x)-ydata));

x0 = [7000,  35, .0001, .25, 0, 0, 0];
lb = [100, 0, .000000007, .25, -1, -1, -0];
ub = [70000, 500, .7, 1, +1, +1, +0];

options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display', 'iter', 'MaxFunEvals', 10000, 'MaxIter', 10^3 );
%options.FunctionTolerance = 1e-10;
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(Yfun, x0, xdata, ydata, lb, ub, options);
 
y = infection(x, xdata, clim_data);
plot([y(2:end) ydata(2:end)]);
yvals = y(2:end);
xvals = ydata(2:end);

%beta = a(3) * (1 - a(5) * precip - a(6) * temp);
%plot(beta.*y(:,1)'/x(4))
%plot(residual2)

rmse = (sum((y(2:end) - ydata(2:end)).^2)/numel(y(2:end))) ^ 0.5;
fprintf('\n\nRMSE = %.3f\n\n',rmse);
RE = (y(2:end) - ydata(2:end))./y(2:end);

%CI = nlparci(x,residual,'jacobian',jacobian, 'alpha', .05);
%J= jacobian;
%  xCovariance = (J.'*J)\residual;

%AIC
RSS = sum((y(2:end) - ydata(2:end)).^2);
AIC = numel(ydata) * log(RSS / numel(ydata)) + 12;
  
[maxval1,ind1] = max(y(2:end));
[maxval2,ind2] = max(ydata(2:end));
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  fid = fopen('Result/data_14.dat', 'w');
  for i=1:numel(xvals)
      fprintf(fid, '%f %f %f \n',10+i,  xvals(i), yvals(i) );
  end
fclose(fid);
  
  fid1 = fopen('Result/param_14.dat', 'w');
  
 fprintf(fid1, 'x0 = ( %f, %f, %f, %f, %f, %f, %f, %f)\n\n', x0(1), x0(2), x0(3), x0(4), x0(5), x0(6), x0(7), x0(8));
 
  fprintf(fid1, '-----------------------------------------------------------------------------\n\n');
  
 fprintf(fid1, 'x = ( %f, %f, %f, %f, %f, %f, %f, %f)\n\n', x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8));


  fprintf(fid1, 'RMSE = %f \n\n', rmse);
  
  fprintf(fid1, 'AIC = %f \n\n',  AIC);
  
  fprintf(fid1, '(Pick Time, Pick Magnitude) = (%f,%f) \n\n', ind1, maxval1);
  
  fprintf(fid1, '(Original Pick Time, Original Pick Magnitude) = (%f, %f) \n\n', ind2, maxval2);

 
fclose(fid1);  
  
%%% G
%%% Genetic Algorithm
%opts = optimoptions('ga','PlotFcn',@gaplotbestf);
%Call the ga solver where x(1) has integer values
%rng(1,'twister') % for reproducibility
%IntCon = 1;
%[x,fval,exitflag] = ga(fun,length(lb),[],[],[],[],lb,ub,[],IntCon,opts);



%  fid2 = fopen('result2/data_35.dat', 'w');
%  for i=1:numel(t)
%      fprintf(fid2, '%f %f %f %f\n',10 + t(i),  y(i,1), y(i,2), y(i,3));
%  end
%fclose(fid2);

%
%   fid = fopen('result1/RMSE2.dat', 'w');
%for i=1:numel(xvals)
%      fprintf(fid, '%f %f\n',i+10, RE(i));
%       end
%  fclose(fid);
%}
