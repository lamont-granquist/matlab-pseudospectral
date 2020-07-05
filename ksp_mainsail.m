close all; clear classes; clc;
format shortG;
format compact;

K = 1; % h steps
M = 40; % p steps

body = Body.Kerbin();

launcher = Launch(body);

initialMass = 43500;

scales = launcher.InitScales(initialMass);

%stages(1) = Stage.StageFromKSP2(scales, 7.00, initialMass, 11500, 310.55, 284.96); % 3105.8
%stages(1) = Stage.StageFromKSP2(scales, 6.00, initialMass, 11500, 310.55, 284.96); % 3106.5
stages(1) = Stage.StageFromKSP2(scales, 5.00, initialMass, 11500, 310.55, 284.96); % 3102.2
%stages(1) = Stage.StageFromKSP2(scales, 4.50, initialMass, 11500, 310.55, 284.96); % 3106.1
%stages(1) = Stage.StageFromKSP2(scales, 4.00, initialMass, 11500, 310.55, 284.96); % 3104.8
%stages(1) = Stage.StageFromKSP2(scales, 3.75, initialMass, 11500, 310.55, 284.96); % 3106.2
%stages(1) = Stage.StageFromKSP2(scales, 3.00, initialMass, 11500, 310.55, 284.96); % 3119.6
%stages(1) = Stage.StageFromKSP2(scales, 2.00, initialMass, 11500, 310.55, 284.96); % 3201
%stages(1) = Stage.StageFromKSP2(scales, 1.90, initialMass, 11500, 310.55, 284.96); % 3216.2
%stages(1) = Stage.StageFromKSP2(scales, 1.80, initialMass, 11500, 310.55, 284.96); % 3235.8
%stages(1) = Stage.StageFromKSP2(scales, 1.70, initialMass, 11500, 310.55, 284.96); % 3261.3
%stages(1) = Stage.StageFromKSP2(scales, 1.60, initialMass, 11500, 310.55, 284.96); % 3295
%stages(1) = Stage.StageFromKSP2(scales, 1.50, initialMass, 11500, 310.55, 284.96); % 3339.4
%stages(1) = Stage.StageFromKSP2(scales, 3.52, initialMass, 11500, 310.55, 284.96); % 3109.8
stages(1).coasts = 1;
stages(1).free_final_time = true;
stages(1).overburn = true;

stages = launcher.FixStages(stages);

auxdata = launcher.InitAuxdata((1.25^2)*pi, 0.5, 5000);

Initial.LatLng(launcher, deg2rad(0), 0);

smaT  = body.radius + 100000;
eccT  = 0.00;
incT  = deg2rad(0);
LANT  = deg2rad(0);
ArgPT = deg2rad(0);
Terminal.Kepler5(launcher, smaT, eccT, incT, LANT, ArgPT);

rT = body.radius + 100000;
vT = sqrt(body.mu/rT);
%Terminal.FPA3(launcher, rT, vT ,0);

launcher.Solve(stages, K, M);
