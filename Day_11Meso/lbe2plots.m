%% Import data from CSV text files

% Set size of lattice used in calculations (needed to work out where each data set starts and ends)

Ymax = 131;

% Import the data from the file and then pull out individual datasets
% (trick is that the file reader adds a line of zeros for any lines without data,
% so we know the starting rows and finishing rows for each set of data)

data = csvread("density.csv", 3, 0);

density = data(1:Ymax,:);
rhoN = data(Ymax+2:2*Ymax+1,:);
kurv = data(2*Ymax+3:3*Ymax+2,:);
psi = data(3*Ymax+4:4*Ymax+3,:);

figure('name','density surface plot');
surf(density);
pause ();

figure('name','density map');
contourf(density,30); axis equal; colorbar;
pause ();

figure('name','phase index surface plot');
surf(rhoN);
pause ();

figure('name','phase index map');
contourf(rhoN,29); axis equal; colorbar;
pause ();

figure('name','curvature surface plot');
surf(kurv);
pause ();

figure('name','curvature map');
contourf(kurv,29); axis equal; colorbar;
pause ();


figure('name','velocity modulus surface plot');
surf(psi);
pause ();

figure('name','velocity modulus map');
contourf(psi,29); axis equal; colorbar;
pause ();



