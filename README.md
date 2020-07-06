# matlab-pseudospectral

Kinda messsy matlab code of me exploring writing a GPOPS-like pseudospectral optimizer in Matlab.

* delta3.m is probably the best entrypoint, that's the classic delta-III launch example from Betts.
* it calls Launch which converts it to the pseudospectral problem.
* then it calls Solver which is the main solver code.
* requires ipopt
* requires casdadi for AD and uses its wrapper around ipopt

Snapshot in time of my messing around exploring stuff and learning math so you'll find exploration dead ends and buggy code.
