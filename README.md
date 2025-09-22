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


# Compilation
This should be compiled using PyInstaller on the same operating system as you intend on running the program. 

# Setting it up
Some template modificaitons are to be done for your own particular run case. The setup was built to follow a certain folder structure as shown here:
```bash
... PARENT-DIR
    |
    |---

```

