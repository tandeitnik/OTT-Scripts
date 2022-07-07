%this scripts evaluates the optical forces along the x,y, and z direction of a trapped spherical particle on an optical tweezers. It also calculates
%the potential well along the x direction. All the forces and the potential well are plotted at the end of the script.

% Add the toolbox to the path (assuming we are in ott/examples)
addpath('../');

% Close open figures
close all;

%% Describe the particle, beam and surrounding medium

% Make warnings less obtrusive
ott.warning('once');
ott.change_warnings('off');

% Refractive index of particle and medium
n_medium = 1.00;
n_particle = 1.46;

% Wavelength of light in vacuum [m]
wavelength0 = 780e-9;

% Calculate the wavelength in the medium
wavelength_medium = wavelength0/n_medium;

% Radius of particle
radius = 143e-9/2;

% Specify the numerical aparture of the beam/microscope
NA = 0.68;

c = 3e8; %velocity of light in vacuum [m/s]
Kb = 1.380649e-23; %Boltzmann constant SI
T = 293; %bath temperature [K]


% Specify the OT power [W]
P = 1;
force_factor = (n_medium*P)/c;

% Create a T-matrix for a sphere
T_matrix = ott.Tmatrix.simple('sphere', radius, 'wavelength0', wavelength0, ...
    'index_medium', n_medium, 'index_particle', n_particle);


% Create a simple Gaussian beam
beam = ott.BscPmGauss('NA', NA, 'polarisation', [ 1 0 ], ...
    'index_medium', n_medium, 'wavelength0', wavelength0);

z = [0;0;1]*linspace(-8,8,80)*wavelength_medium;
fz = ott.forcetorque(beam, T_matrix, 'position', z);

zeq = ott.find_equilibrium(z(3, :), fz(3, :));
if isempty(zeq)
  warning('No axial equilibrium in range!')
  zeq=0;
end
zeq = zeq(1);

% Calculate force along x-axis (with z = zeq, if found)
r = [1;1;0]*linspace(-8,8,500)*wavelength_medium + [0;0;zeq];
fr = ott.forcetorque(beam, T_matrix, 'position', r);

x = r(1,:);
y = r(2,:);
z = z(3,:);
fx = fr(1,:)*force_factor;
fy = fr(2,:)*force_factor;
fz = fz(3,:)*force_factor;
potential_x = (-1*cumtrapz(x,fx))/(Kb*T);

%plotting some relevant graphs

%the potential along x
figure(1); plot(x/wavelength_medium,potential_x);
xlabel('{\it x} [\lambda_m]');
ylabel('{\it E} [K_B T]');
aa = axis;
hold on;
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
set(gca,'FontSize',25)
hold off;
%the force along z
figure(2); plot(z/wavelength_medium,fz*1e12);
xlabel('{\it z} [\lambda_m]');
ylabel('{\it F_z} [pN]');
aa = axis;
hold on;
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
set(gca,'FontSize',25)
hold off;
%the force along x
figure(3); plot(x/wavelength_medium,fx*1e12);
xlabel('{\it x} [\lambda_m]');
ylabel('{\it F_x} [pN]');
aa = axis;
hold on;
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
set(gca,'FontSize',25)
hold off;
%the force along y
figure(4); plot(y/wavelength_medium,fy*1e12);
xlabel('{\it y} [\lambda_m]');
ylabel('{\it F_y} [pN]');
aa = axis;
hold on;
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
set(gca,'FontSize',25)
hold off;
