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

c = 3e8; %velocity of light in vacuum [m/s]
Kb = 1.380649e-23; %Boltzmann constant SI
T = 293; %bath temperature [K]

% Wavelength of light in vacuum [m]
wavelengthCentral = 780e-9;
freqCentral = c/wavelengthCentral;
deltaFreq = 50e9;
n = 100;
freqArray = linspace(freqCentral - deltaFreq,freqCentral + deltaFreq,n);
deltaArray = linspace(0-deltaFreq,0+deltaFreq,n);

points = 200;

z = [0;0;1]*linspace(-4,4,points)*(wavelengthCentral/n_medium);
r = [1;0;0]*linspace(-4,4,points)*(wavelengthCentral/n_medium);

% Radius of particle
radius = 143e-9/2;

% Specify the numerical aparture of the beam/microscope
NA = 0.77;
    
% Specify the OT power [W]
P = 0.6;

fxArray = zeros(n,points);
fzArray = zeros(n,points);

for i = 1:n
    
    wavelength0 = c/freqArray(i);

    % Calculate the wavelength in the medium
    wavelength_medium = wavelength0/n_medium;
    
    %force_factor = (n_medium*P)/c;
    force_factor = 1;

    % Create a T-matrix for a sphere
    T_matrix = ott.Tmatrix.simple('sphere', radius, 'wavelength0', wavelength0, ...
        'index_medium', n_medium, 'index_particle', n_particle);
    
    
    % Create a simple Gaussian beam
    beam = ott.BscPmGauss('NA', NA, 'polarisation', [ 1 0 ], ...
        'index_medium', n_medium, 'wavelength0', wavelength0);
    
    
    fz = ott.forcetorque(beam, T_matrix, 'position', z);
    
    zeq = ott.find_equilibrium(z(3, :), fz(3, :));
    if isempty(zeq)
      warning('No axial equilibrium in range!')
      zeq=0;
    end
    zeq = zeq(1);
    
    % Calculate force along x-axis (with z = zeq, if found)
    
    fr = ott.forcetorque(beam, T_matrix, 'position', r);
    
    
    fxArray(i,:) = fr(1,:)*force_factor;
    fzArray(i,:) = fz(3,:)*force_factor;
    %potential_x = (-1*cumtrapz(x,fx))/(Kb*T);

end

stdZ = zeros(1,points);
stdX = zeros(1,points);
meanZ = zeros(1,points);
meanX = zeros(1,points);

for i = 1:points

    stdZ(1,i) = std(fzArray(:,i));
    stdX(1,i) = std(fxArray(:,i));
    meanZ(1,i) = mean(fzArray(:,i));
    meanX(1,i) = mean(fxArray(:,i));

end

errorbar(r(1,:)*1e6,meanX,stdX,"o")
xlabel('x [\mum]')
ylabel('\langleF\rangle')
aa = axis;
hold on;
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
set(gca,'FontSize',25)
hold off;

errorbar(z(3,:)*1e6,meanZ,stdZ,"o")
xlabel('z [\mum]')
ylabel('\langleF\rangle')
aa = axis;
hold on;
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
set(gca,'FontSize',25)
hold off;

%heatmap for z

skip = 30;
yvalues = {};

for i = 1:points

    if mod(i,skip) == 1 || i == n

        yvalues{end+1} = string(z(3,i).*1e6);

    else

        yvalues{end+1} = " ";

    end

end

skip = 13;
xvalues = {};

for i = 1:n

    if mod(i,skip) == 1 || i == n

        xvalues{end+1} = string(deltaArray(1,i)./1e9);

    else

        xvalues{end+1} = " ";

    end

end

h = heatmap(transpose(fzArray),'GridVisible','off');
h.YDisplayLabels = yvalues;
h.XDisplayLabels = xvalues;
h.XLabel = '\Deltaf [GHz]';
h.YLabel = 'z [\mum]';
set(gca,'FontSize',25)


%heatmap for x

skip = 21;
yvalues = {};

for i = 1:points

    if mod(i,skip) == 1 || i == n

        yvalues{end+1} = string(r(1,i).*1e6);

    else

        yvalues{end+1} = " ";

    end

end

skip = 13;
xvalues = {};

for i = 1:n

    if mod(i,skip) == 1 || i == n

        xvalues{end+1} = string(deltaArray(1,i)./1e9);

    else

        xvalues{end+1} = " ";

    end

end

h = heatmap(transpose(fxArray),'GridVisible','off');
h.YDisplayLabels = yvalues;
h.XDisplayLabels = xvalues;
h.XLabel = '\Deltaf [GHz]';
h.YLabel = 'x [\mum]';
