%1st part, simulate traces

kb = 1.38e-23; % Boltzmann cte. in J / K
T = 293.0; % temperature in Kelvin
t_max = 1; % [s] max simulation time
dt = 1e-6; % time interval for the simulation [s]
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
    waitbar(M/trials, wb, sprintf('Progress: %d %%', floor(M/trials*100)));
    
end

close(wb)

simulatedSignal = {};

for i = 1:trials

    simulatedSignal{end+1} = simulationPositions{i}(1:round(dt_s/dt):end);

end

%2st part, apply mean squared displacement method

%first estimating k using the equipartition method
%(the estimated kTentative will be used as a start point on the non-linear
%fitting)
k_array = zeros(1,trials);

for i = 1:trials

    k_array(1,i) = kbT/var(simulatedSignal{i});

end

kTentative = mean(k_array);
tRelax = gamma/kTentative;

%evaluating the MSD
divisions = 100;
tauMin = 0.5*tRelax;
tauMax = 6*tRelax; %[s]
tau = linspace(tauMin,tauMax,divisions);
MSDarray = zeros(1,divisions);
MSDmatrix = zeros(trials,divisions);

for M = 1:trials

    for i = 1:divisions
    
        delta = round(tau(i)/dt);
        positions = simulationPositions{M};
        MSD = 0;
        for j = 1:N-delta
            MSD = MSD + (positions(j+delta)-positions(j))^2;
        end
        MSDarray(1,i) = MSDarray(1,i) + MSD/(N-delta);
        MSDmatrix(M,i) = MSD/(N-delta);
    end

end

MSDarray = MSDarray./trials;

%evaluating the variance
MSDvar = zeros(1,divisions);
for j = 1:divisions
    for i = 1:trials
        MSDvar(1,j) = MSDvar(1,j) + (MSDmatrix(i,j) - MSDarray(1,j))^2;
    end
end
MSDvar = MSDvar./trials;
err = sqrt(MSDvar);

%fitting the MSD to the theorical model
% Define Start points, fit-function and fit curve
k0 = [kTentative]; 
fitfun = fittype( @(k,x) (2*kbT/k)*(1-exp((-1*x*k)/gamma)) );
[fitted_curve,gof] = fit(tau(:),MSDarray(:),fitfun,'StartPoint',k0);
coeffvals = coeffvalues(fitted_curve);
bounds = confint(fitted_curve);

k = coeffvals;
k_error = abs(bounds(1,1) - bounds(2,1));
f_c = k/(2*pi*gamma);
f_c_err = k_error/(2*pi*gamma);

%3rd part, plotting the result

errorbar(tau*1000,MSDarray*1e18,err*1e18,'o','DisplayName','experimental MSD')
hold on
plot(tau*1000,fitted_curve(tau)*1e18,'--r','LineWidth',2,'DisplayName','Non-linear fitting')
xlabel('{\it \tau} [ms]');
ylabel('{\it MSD} [nm^2]');
aa = axis;
set(gca,'FontSize',25)
grid on
xlim([tauMin*1000 tauMax*1000])
ylim([MSDarray(1)*1e18,MSDarray(end)*1e18*1.05])
legend('Location','southeast','Orientation','horizontal')
hold off;