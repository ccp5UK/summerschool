%% Import data from CSV text files

% Import the data
datadensity = csvread("data_density.csv");
datastfn = csvread("data_stfn.csv");

figure('name','density surface plot');
surf(datadensity);
pause ();

figure('name','density contour plot');
contourf(datadensity,30);
pause ();

figure('name','stream function surface plot');
surf(datastfn);
pause ();

figure('name','streamlines (stream function contour plot)');
contourf(datastfn,30);
pause ();
