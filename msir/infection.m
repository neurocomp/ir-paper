% infection.m 
 
function Y = infection(a, xdata, clim_data) 

global itter
fprintf('Itter = %d\n', itter);
fprintf('Solution = (%d, %d, %d, %d) \n', a(1), a(2), a(3), a(4));
fprintf(' ... \n');
itter=itter+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ft = [1:.01:xdata(end)];
  shum = interp1(1: max(xdata), clim_data(:,1), ft);
  temp = interp1(1: max(xdata), clim_data(:,2), ft);
  precip = interp1(1: max(xdata), clim_data(:,3), ft);
  
  ic = [a(1), a(2)];
  opts = odeset('RelTol', 1e-2, 'AbsTol', 1e-4);
  [t,y] = ode45(@(t,x) sir(t, a, x, ft, precip, temp, shum), ft, ic, opts);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Y = zeros(max(xdata), 1);
  for i = (2):(max(xdata))
      %[t,tspan]
      alpha = find(t <= i - 1);
      beta = find(t < i);
      Dt = diff(t);
%      dt = Dt(1);
      %Y(i) = a(5) * dt * sum(y((alpha(end):beta(end)),2));
      Y(i) =a(3) * sum(y(alpha(end):beta(end),1) .*  y(alpha(end):beta(end),2))*.01;
  end
  %U = interp1(t, y(:,2), x);
  
 
  
  
  