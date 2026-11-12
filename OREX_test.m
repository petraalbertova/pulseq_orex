%% Initialization
clear

% loading pulseq objects from .mat file
% contains AHP pulses and spiral readout
% objects calculated via OpenMRF: openmrf.org
load('pulseq_objects.mat');

% define system definitions
system = mr.opts('maxGrad', 35,  'gradUnit', 'mT/m', ...
                 'maxSlew', 170, 'slewUnit', 'T/m/s', ...
                 'B0',      3, ...
                 'rfDeadTime',          1e-4, ...                 
                 'rfRingdownTime',      1e-4, ...
                 'rfRasterTime',        1e-6, ...
                 'adcDeadTime',         1e-5, ...
                 'adcRasterTime',       1e-7, ...
                 'gradRasterTime',      1e-5, ...
                 'blockDurationRaster', 1e-5 );

% init seq object
seq = mr.Sequence(system);

% FOV definitions: slice
fov_dz       =  5 *1e-3; % slice thickness
fov_z_offset = 10 *1e-3; % slice offset; we need a small slice offset for the TREX gradients 

% FOV definitions: in plane
% do not change!
% spiral gradients were generated for the following parameters
Nxy    = 96;
fov_xy = 256 *1e-3; 

% set .seq file definitions
seq.setDefinition('Name', 'OREX_test');
seq.setDefinition('FOV', [fov_xy fov_xy fov_dz]);

%% Spin-Lock parameters & objects
tSL = 100 *1e-3;
tSL = round(tSL / (2*system.rfRasterTime)) * 2*system.rfRasterTime;
fSL = 120;

% spin-lock object
SL        = mr.makeBlockPulse(pi, system, 'Duration', tSL/2, 'freqOffset', 0, 'PhaseOffset', 0, 'use', 'preparation'); % contains dummy flipangle pi
SL.signal = SL.signal*0 + fSL;

% adding small delays between consecutive RF pulses
d10us = mr.makeDelay(1e-5);

%% Slice Selective Refocusing

% refocusing pulse
% note: for the paper, we used a sigpy slr pulse with the following settings:
% alpha=pi, duration=5e-3, tbw=3, ptype='se', ftype='ls', d1=0.01, d2=0.01, cancel_alpha_phs=false
[rf, gz] = mr.makeSincPulse( pi, ...
                             system, ...
                             'Duration', 5e-3, ...
                             'timeBwProduct', 3, ...
                             'apodization', 0.5, ...
                             'SliceThickness', fov_dz, ...
                             'use', 'refocusing');
rf.freqOffset = fov_z_offset * gz.amplitude; % apply frequency offset for slice shift -> TREX gradients

% crusher parameters
ntwist_xy = 8;
ntwist_z  = 8;

% calculate crusher gradients
gx_crush = mr.makeTrapezoid('x', 'Area', (ntwist_xy/2) / (fov_xy/Nxy), 'system', system);
gy_crush = mr.makeTrapezoid('y', 'Area', (ntwist_xy/2) / (fov_xy/Nxy), 'system', system);
gz_crush = mr.makeTrapezoid('z', 'Area', (ntwist_z/2)  / fov_dz - gz.area/2, 'Duration', mr.calcDuration(gx_crush), 'system', system);

%% TREX parameters
% the TREX gradients emulate the biomagnetic field in the offcenter slice

% magnetic field oscillation parameters
dB0      = 10 *1e-9; % [T] emulate a 10nT field oscillation in the 10mm offcenter slice
gtrex    = dB0 * system.gamma / fov_z_offset; % [Hz/m] corresponding gradient amplitude
freq     = 120; % [Hz] oscillation frequency; resonant condition if fSL=freq
n_phase  = 20; % number of relative phases between spin-lock and field oscillation
phi      = linspace(0, 2*pi, n_phase+1);
phi(end) = [];

% recovery time
trec = mr.makeDelay(1);

%% Create sequence

% add 3 dummy scans for steady-state
n_dummy = 3;

for loop_phase = 1-n_dummy : n_phase

    % calculate TREX gradients for current phase
    temp_idx = max(1, loop_phase);
    temp_t   = (0:system.gradRasterTime:tSL-system.gradRasterTime)';
    temp_g   = gtrex * sin(2*pi*freq*temp_t + phi(temp_idx));    
    temp_n   = numel(temp_g)/2;
    gtrex_1  = mr.makeArbitraryGrad('z', temp_g(1:temp_n),     'first', temp_g(1),        'last', temp_g(temp_n), 'system', system);
    gtrex_2  = mr.makeArbitraryGrad('z', temp_g(temp_n+1:end), 'first', temp_g(temp_n+1), 'last', temp_g(end),    'system', system);    

    % adiabatic tip down pulse (+z -> +x)
    seq.addBlock(ahp1);
    seq.addBlock(d10us);

    % spin-lock part 1
    seq.addBlock(SL, gtrex_1);
    seq.addBlock(d10us);

    % spin-lock part 2
    seq.addBlock(SL, gtrex_2);
    seq.addBlock(d10us);

    % adiabatic tip down pulse (+x -> -z)
    seq.addBlock(ahp2);

    % slice selective refocusing and crusher gradients
    seq.addBlock(gx_crush, gy_crush, gz_crush);
    seq.addBlock(rf, gz);
    seq.addBlock(gx_crush, gy_crush, gz_crush);

    % spiral readout
    if loop_phase < 1
        seq.addBlock(gx_read, gy_read, gz_spoil);
    else
        seq.addBlock(gx_read, gy_read, gz_spoil, adc);
    end

    % recovery delay
    seq.addBlock(trec);

end

%% check timings, plot sequence, export .seq file
[timings_ok, error_report] = seq.checkTiming;
seq.plot('stacked', true, 'showGuides', false);
seq.write('OREX_test.seq', true);

%% calculate k-space trajectory for reconstruction

% init clean seq object
seq_traj = mr.Sequence(system);

% write dummy seqeunce
seq_traj.addBlock(mr.makeSincPulse( pi/2, system, 'Duration', 2e-3, 'timeBwProduct', 4, 'use', 'excitation'));
seq_traj.addBlock(gx_read, gy_read, adc);

% calculate k-space
ktraj = seq_traj.calculateKspacePP();
