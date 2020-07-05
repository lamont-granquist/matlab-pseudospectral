close all; clear classes; clc;
format shortG;
format compact;

K = 1; % h steps
M = 40; % p steps

body = Body.Earth();

launcher = Launch(body);

initialMass = 58171;

scales = launcher.InitScales(initialMass);

stages(1) = Stage.StageFromKSP(scales, 1.46, 1.30, 34227, 3974, 150.1);
stages(2) = Stage.StageFromKSP(scales, 1.64, 1.41, 1322, 807, 54.9);
stages(3) = Stage.StageFromKSP(scales, 3.30, 2.83, 657, 142, 54.9);
%stages(1) = Stage.StageFromKSP(scales, 1.65, 1.30, 58171, 3979, 172.9);
%stages(2) = Stage.StageFromKSP(scales, 8.61, 7.62, 400, 241, 12.5);
%stages(1) = Stage.StageFromKSP(scales, 1.60, 1.26, 59905, 5664, 172.9);
%stages(2) = Stage.StageFromKSP(scales, 1.62, 1.43, 2127, 791, 105.0);
%stages(3) = Stage.StageFromKSP(scales, 4.62, 3.96, 470, 114, 37.8);
stages(end).free_final_time = true;
stages(end).overburn = true;
stages(end).coasts = 1;
%stages(end).coast_before = true;

%bt = 150.1         54.9       51.533

stages = launcher.FixStages(stages);

auxdata = launcher.InitAuxdata((1.25^2)*pi, 0.5, 5000);

Initial.LatLng(launcher, deg2rad(28.5), 0);

PeA = 185000;
ApA = 185000;

PeR = PeA + body.radius;
ApR = ApA + body.radius;

smaT = ( PeR + ApR ) / 2;
eccT = ( ApR - PeR ) / ( ApR + PeR );
incT  = deg2rad(28.5);
LANT  = deg2rad(270);
ArgPT = deg2rad(97.686);
Terminal.Kepler5(launcher, smaT, eccT, incT, LANT, ArgPT);

%launcher.Solve(stages, K, M);

rT = PeR;
vT = sqrt(body.mu * (2/rT - 1/smaT));
%Terminal.FPA3(launcher, rT, vT, 0);
%Terminal.FPA4(launcher, rT, vT, 0, incT);
%Terminal.FPA5(launcher, rT, vT, 0, incT, LANT);

launcher.Solve(stages, K, M);
