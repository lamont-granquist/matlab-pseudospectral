close all; clear classes; clc;
format longG;
format compact;

K = 1; % h steps
M = 20; % p steps

body = Body.Earth();

launcher = Launch(body);

srbWetMass = 19290;
srbPropMass = 17010;
srbDryMass = srbWetMass - srbPropMass;
srbBurnTime = 75.2;
srbThrust = 628500;
srbMdot = srbPropMass / srbBurnTime;

firstWetMass = 104380;
firstPropMass = 95550;
firstDryMass = firstWetMass - firstPropMass;
firstBurnTime = 261;
firstThrust = 1083100;
firstMdot = firstPropMass / firstBurnTime;

secondWetMass = 19300;
secondPropMass = 16820;
secondDryMass = secondWetMass - secondPropMass;
secondBurnTime = 700;
secondThrust = 110094;
secondMdot = secondPropMass / secondBurnTime;

payloadMass = 4164;

initialMass = payloadMass + secondWetMass + firstWetMass + 9 * srbWetMass;

scales = launcher.InitScales(initialMass);

stages(1) = Stage.StageFromMdot(scales, 6 * srbThrust + firstThrust, initialMass, 6 * srbMdot + firstMdot, srbBurnTime);
stages(1).SetAe( exit_area(firstMdot, 302, 255) + 6 * exit_area(srbMdot, 279.8, 251) );

stages(2) = Stage.StageFromMdot(scales, 3 * srbThrust + firstThrust, (stages(1).mf - (6 * srbDryMass / scales.mass)) * scales.mass, 3 * srbMdot + firstMdot, srbBurnTime);
stages(2).SetAe( exit_area(firstMdot, 302, 255) + 3 * exit_area(srbMdot, 279.8, 251) );

stages(3) = Stage.StageFromMdot(scales, firstThrust, (stages(2).mf - (3 * srbDryMass / scales.mass)) * scales.mass, firstMdot, firstBurnTime - 2 * srbBurnTime);
stages(3).SetISP(302, 255);

stages(4) = Stage.StageFromMdot(scales, secondThrust, (stages(3).mf - firstDryMass / scales.mass ) * scales.mass, secondMdot, secondBurnTime);
stages(4).SetAe( 0 );
stages(4).free_final_time = true;
%stages(4).coast_before = true;
%stages(4).coasts = 1;

stages = launcher.FixStages(stages);

auxdata = launcher.InitAuxdata(4*pi, 0.5, 2000);

Initial.LatLng(launcher, deg2rad(28.5), 0);

smaT  = 24361140;
eccT  = 0.7308;
incT  = deg2rad(28.5);
LANT  = deg2rad(269.8);
ArgPT = deg2rad(130.5);
Terminal.Kepler5(launcher, smaT, eccT, incT, LANT, ArgPT);

rT = body.radius + 1000;
vT = sqrt(body.mu/rT);
%Terminal.FPA3(launcher, rT, vT ,0);

launcher.Solve(stages, K, M);

