% this scripts simulates the dynamics of a particle trapped on an optical tweezers and uses the potential analysis method to evaluate the
% spring stifiness k of the trap for the x axis. Can be easily modified to evaluate for the other axis

kb = 1.38e-23; % Boltzmann cte. in J / K
T = 293.0; % temperature in Kelvin
t_max = 1; % [s] max simulation time
dt = 1e-5; % time interval for the simulation [s]
dt_s = 1e-3; %sampling dt [s]
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
trials = 50;

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

%determining x_max and x_min of all traces
x_max = max(simulationPositions{1});
x_min = min(simulationPositions{1});

for i = 2:trials

    x_max = max(x_max,max(simulationPositions{i}));
    x_min = min(x_min,min(simulationPositions{i}));

end

%building the edges
nbins = 50;
edges = linspace(x_min,x_max,nbins);
edgesStamps = edges(1:end-1);

meanU = -kbT.*log(histcounts(simulationPositions{1},edges, 'Normalization', 'probability'));

for i = 2:trials

    meanU = meanU + -kbT.*log(histcounts(simulationPositions{i},edges, 'Normalization', 'probability'));

end

meanU = meanU/trials;

varU = (-kbT.*log(histcounts(simulationPositions{1},edges, 'Normalization', 'probability')) - meanU).^2;

for i = 2:trials

    varU = varU + (-kbT.*log(histcounts(simulationPositions{i},edges, 'Normalization', 'probability')) - meanU).^2;

end

varU = varU/(trials-1);

%fitting meanU to a second degree polynomial
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( edgesStamps, meanU );

% Set up fittype and options.
ft = fittype( 'poly2' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
bounds = confint(fitresult);

k = 2*fitresult.p1;
k_error = (bounds(1,1) - bounds(1,2))/2;
f_c = k/(2*pi*gamma);
f_c_err = k_error/(2*pi*gamma);
