% this scripts simulates the dynamics of a particle trapped on an optical tweezers and uses the equipartition method to evaluate the
% spring stifiness k of the trap for the x axis. Can be easily modified to evaluate for the other axis

kb = 1.38e-23; % Boltzmann cte. in J / K
T = 293.0; % temperature in Kelvin
t_max = 1; % [s] max simulation time
dt = 1e-5; % time interval for the simulation [s]
dt_s = 1e-3; %sampling dt [s] - should be greater then the relaxation time of the trap so the data points are uncorrelated
radius = 1.15e-6/2; % Radius of particle [m]
viscosity = 0.0008538; % water [N/m^2] https://www.omnicalculator.com/physics/water-viscosity
%viscosity = 1.6e-5; %air
gamma =  6*pi*viscosity*radius; % # damping coef. [N.s/m]
n_medium = 1.33; % Medium refractive index
n_particle = 1.46; % Particle refractive index
wavelength0 = 780e-9; % Wavelength of light in vacuum [m]
c = 299792458; %Speed of light [m/s]
P = 30e-3; %Laser power [W]
NA = 1.3; %Numerical Aperture
trials = 50; %number of traces that will be generated, i.e., number of experiments

kbT = kb*T;
wavelength_medium = wavelength0/n_medium;

% Create a T-matrix for a sphere
T_matrix = ott.Tmatrix.simple('sphere', radius, 'wavelength0', wavelength0, ...
    'index_medium', n_medium, 'index_particle', n_particle);
% Create a simple Gaussian beam
beam = ott.BscPmGauss('NA', NA, 'polarisation', [ 1 0 ], ...
        'index_medium', n_medium, 'wavelength0', wavelength0);

force_factor = n_medium*P/c;

timestamps = linspace(0,t_max,t_max/dt);
N = size(timestamps,2);
simulationPositions = {};
simulatedSignal = {};

z = [0;0;1]*linspace(-10,10,500)*wavelength_medium;
fz = ott.forcetorque(beam, T_matrix, 'position', z);

% Find the equilibrium along the z axis
zeq = ott.find_equilibrium(z(3, :), fz(3, :));

r = [1;1;0]*linspace(-10,10,500)*wavelength_medium + [0;0;zeq];
fr = ott.forcetorque(beam, T_matrix, 'position', r);

x = r(1,:);
y = r(2,:);
fx = fr(1,:);
fy = fr(2,:);

%to save time, instead of using OTT to evaluate the force at each time step, I consider that the problem 1D by ignoring the y and z movements (considering
%that the particle keeps at y = 0 and z = zeq) and by fitting the force along x by a smooth spline. I use the fitted sliple to evaluate the force at each 
#time step, which is some order of magnitudes faster then calling OTT's function. The result is consistent!

%% Fit: 'force_x'.
[xData, yData] = prepareCurveData( x, fx );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );

% Fit model to data.
[fit_fx, gof] = fit( xData, yData, ft );

wb = waitbar(0, 'Starting');

for M = 1:trials

    positions = zeros([1,N]);
    
    for i = 2:N
       
        f = fit_fx(positions(1,i-1))*force_factor;
        W = sqrt(2.0 * kbT * dt / gamma) * normrnd(0,1,[1,1]);
        positions(1,i) = positions(1,i-1) +f*dt/gamma + W;
       
    
    end

    simulationPositions{end+1} = positions;
    simulatedSignal{end+1} = positions(1:dt_s/dt:end);
    waitbar(M/trials, wb, sprintf('Progress: %d %%', floor(M/trials*100)));
    
end

close(wb)

k_array = zeros(1,trials);

for i = 1:trials

    k_array(1,i) = kbT/var(simulatedSignal{i});

end

k = mean(k_array); #trap stifiness along x
k_err = std(k_array); #error of k
f_c = k/(2*pi*gamma); #corner frequency if a PSD where to be taken
f_c_err = k_err/(2*pi*gamma); #error of f_c

