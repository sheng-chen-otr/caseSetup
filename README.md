# INTRODUCTION
This is a program to setup an OpenFOAM case with ease. It is mainly built around external aerodynamics for automotive applications, however it is capable of some internal flow applications as well. 

# Capabilities

## Compatabilities
OpenFOAM version compabilities: v2012 - v2212 (newer versions untested but should work)

Operating system compatilibies: Ubuntu 22.XX, might work on Windows but untested

## Solvers

| Solvers/Utility  | Turbulence Modeling Method | Turbulence Models |
| -------------    | -------------              | -------------     |
| SimpleFoam       | RANS                       | KO-SST, SA        |
| PisoFoam         | DDES                       | KO-SST, SA        |
