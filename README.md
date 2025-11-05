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
    |--- 01_INCOMING_DATA
    |--- 02_REFERENCES
    |       |
    |       |--- MSH #this is the location of your trisurfaces
    |       |--- CAD
    |
    |--- CASES
            |
            |--- 001
            |--- 002
            |--- 003 #these are you cases
                  |--- caseSetup
```



# Quick Start
After compiling to an executable, the templates are located inside '_internal'. The case templates are located in '_internal/setupTemplates'. 

Executing 'caseSetup' with no flags will result in the use of the default template. This template is setup for external aerodynamics. Additional templates can be added or made by copying the 'default' template and making your modifications within.

This version, v4.x, is not compatible with previous versions, as there was a major overhaul in the process.

If you have never ran this version or later you would have to have 'caseSetup' generate a new 'caseSetup' file. Fill in the 'caseSetup' file with the appropriate settings. If using all default settings, just fill in the geometry section. 

This is an example of the geometry. Each geometry makes up one line. The order is as such 'geometry.obj.gz, SCALE, REFINEMENT LEVEL, NLAYERS, RELATIVE LAYER GROWTH, WALL MODEL'

'WALL MODEL [OPTIONS: high, low]' refers to high or low Reynolds number model. High Reynolds number would mean the flow velocity is very high and a wall model is used, low Reynolds would not use any wall modell and resolve the boundary layer directly. 


```bash
[GEOMETRY]
GEOM = FRONT_WING.obj.gz,1,8,5,1.1,high
        REAR_WING.obj.gz,1,8,5,1.1,high
        BODY.obj.gz,1,8,5,1.1,high
        ROTA-FRONT-WHEEL_LHS.obj.gz,1,8,5,1.1,high
        ROTA-FRONT-WHEEL_RHS.obj.gz,1,8,5,1.1,high
        ROTA-REAR-WHEEL_LHS.obj.gz,1,8,5,1.1,high
        ROTA-REAR-WHEEL_RHS.obj.gz,1,8,5,1.1,high



```
