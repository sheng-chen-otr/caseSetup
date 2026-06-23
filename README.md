# caseSetup-dev: OpenFOAM Case Setup & Suspension Kinematics Framework

## Overview

**caseSetup** is a Python-based framework for automated OpenFOAM case setup with advanced support for external aerodynamics and suspension kinematics. It streamlines mesh generation, boundary condition configuration, solver setup, and post-processing workflows for automotive CFD applications.

**Key Feature**: **Double Wishbone + Pushrod Suspension Kinematic Solver** — automatically computes per-suspension-component transforms (translation, rotation, pivot point) from ride-height inputs, enabling accurate geometry deformation that follows actual suspension mechanics.

---

## Compatibility

| Aspect | Support |
|--------|---------|
| **OpenFOAM Versions** | v2012 - v2212 (newer versions untested) |
| **Operating Systems** | Ubuntu 22.XX (Linux preferred); Windows/macOS untested |
| **Python** | 3.8+ |

---

## Directory Structure

```
caseSetup-dev/
├── caseSetup.py                    # Main entry point (command-line interface)
├── utilities.py                    # Core geometry/suspension utilities
├── writeSystem.py                  # OpenFOAM system/ directory generators
├── writeConstant.py                # OpenFOAM constant/ directory generators
├── writeScripts.py                 # Solver and cluster script generators
├── testing.py                      # Development/testing utilities
├── compileFiles.sh                 # PyInstaller compilation script
├── README.md                       # This file
│
├── setupTemplates/                 # Case setup templates
│   ├── default/                    # Standard automotive aerodynamics setup
│   │   ├── defaultANSA/            # ANSA mesh configuration
│   │   ├── defaultBCTemplates/     # Boundary condition templates
│   │   ├── defaultCluster/         # Cluster/HPC submission configs
│   │   ├── defaultDicts/           # OpenFOAM dict templates
│   │   ├── defaultPost/            # Post-processing setup
│   │   ├── defaultRefinements/     # Mesh refinement definitions
│   │   └── defaultSetup/           # Ride-height & general config
│   └── otr/                        # Optimal Test Run setup variant
│       ├── defaultANSA/
│       ├── defaultBCTemplates/
│       ├── defaultCluster/
│       ├── defaultDicts/
│       ├── defaultPost/
│       ├── defaultRefinements/
│       ├── defaultFidelity/        # Fidelity-specific settings
│       └── defaultSetup/
│
├── zeroTemplates/                  # OpenFOAM time-zero (initial condition) templates
│   ├── boundaryConditions/         # BC definitions for walls, inlets, outlets
│   ├── models/                     # Turbulence model templates
│   └── solvers/                    # Solver-specific templates (simpleFoam, pimpleFoam, etc.)
│
├── meshUtilities/                  # Mesh generation & analysis tools
│   ├── ansaMesh.py                 # ANSA mesh import/export utilities
│   ├── fidelityMesh.py             # Fidelity mesh integration
│   ├── getStats.py                 # Mesh quality statistics
│   └── _obsolete/                  # Legacy mesh scripts
│
├── postUtilities/                  # Post-processing & results analysis
│   ├── pvPost.py                   # ParaView automation & field extraction
│   ├── postRun.py                  # Post-run analysis orchestration
│   ├── summary.py                  # Results summary generation
│   ├── plotForces.py               # Force coefficient visualization
│   ├── forceConvergencePlot.py     # Convergence analysis
│   ├── createMovies.py             # Animation generation
│   ├── pptGeneration.py            # PowerPoint report generation
│   ├── resultLogWriter.py          # Results logging
│   ├── gsheetExport.py             # Google Sheets export
│   ├── estimateStatisticalError.py # Uncertainty quantification
│   ├── changeCoeff.py              # Force coefficient post-processing
│   ├── postProReportGen.py         # Report generation orchestrator
│   ├── postProReportTemplate/      # LaTeX report templates
│   ├── pptGeneration/              # PowerPoint templates
│   ├── default/                    # Default config & setup files
│   └── _obsolete/                  # Legacy post-processing scripts
│
├── preUtilities/                   # Pre-processing utilities
│   ├── archiveCase.py              # Case archival/backup utilities
│   └── copyTrial.py                # Trial case copying/cloning
│
├── rideHeightUtils/                # **Suspension Kinematics System**
│   ├── doubleWishbonePushrodKinematics.py  # Core kinematic solver
│   ├── suspensionHardpoints_template.cfg   # Hardpoint geometry template
│   ├── suspensionHardpoints_template.json  # JSON format alternative
│   ├── rideHeightConfig                    # Ride-height configuration
│   └── rideHeightMorph.py                  # Suspension morphing utilities
│
├── documentation/                  # Project documentation
│   └── source/                     # LaTeX/technical documentation
│
└── addTemplates/                   # (Not shown) Additional custom templates
```

---

## Core Modules

### **caseSetup.py** — Main Entry Point

Command-line interface for case setup. Orchestrates geometry linkage, mesh configuration, solver setup, and optional ride-height processing.

**Usage**:
```bash
# Standard case setup
python caseSetup.py -s otr

# Create new caseSetup template
python caseSetup.py --new

# Ride-height study: set RUN_RIDE_HEIGHT = True in caseSetup, then run normally.
# caseSetup builds the child cases and re-invokes itself with --rideHeightMode
# internally (do not pass --rideHeightMode yourself).
python caseSetup.py -s otr

# Options:
#  -s, --setup            Setup template type (default: 'default')
#  -d, --controlDict      Write controlDict only
#  --new                  Create new caseSetup file
#  --modules              Show available modules
#  --postProDict          Copy post-processing config
#  --rideHeightMode       INTERNAL ONLY (used by caseSetup when re-running child cases)
```

### **utilities.py** — Core Geometry & Kinematics Engine

**Key Responsibilities**:
- Suspension hardpoint loading & parsing (CFG/INI format)
- Kinematic solver integration (DoubleWishbonePushrodSolver)
- **Per-component transform computation** (WHEEL, UCA, LCA, ROCKER, PUSHROD, DAMPER, TIE)
- Geometry transformation (OBJ/STL with PID preservation)
- Ride-height calculation and case generation

**Key Functions**:
- `calculateRideHeights()` — Compute pitch, roll, heave, wheel movements from ride-height input file
- `_load_suspension_kinematics_setup()` — Load hardpoints and validate kinematic setup
- `_compute_component_transforms()` — **[NEW]** Compute per-component translation/rotation/pivot from kinematic solver state
- `_categorize_pid_keywords_by_component()` — **[NEW]** Match geometry PIDs to suspension components
- `transformGeom()` — Apply per-component transforms to geometry files
- `transformGeometryByPIDRegex()` — PID-aware geometry transformation with regex matching

### **writeSystem.py** — OpenFOAM System Directory

Generates `system/` directory contents:
- `blockMeshDict` — Structured mesh parameters
- `snappyHexMeshDict` — Mesh refinement & snapping
- `decomposeParDict` — Parallel decomposition
- `controlDict` — Solver time stepping & output
- `fvSchemes`, `fvSolution` — Discretization & solver schemes

### **writeConstant.py** — OpenFOAM Constant Directory

Generates `constant/` directory:
- Turbulence model setup
- Material properties
- Geometry boundary condition patches
- Transport properties

### **writeScripts.py** — Job Submission & Execution Scripts

Generates cluster submission scripts:
- SLURM, PBS, SGE configurations
- Meshing, solving, post-processing pipelines
- Job monitoring and result export

### **meshUtilities/** — Mesh Tools

- **ansaMesh.py**: ANSA Hypermesh integration for mesh import/export
- **fidelityMesh.py**: Mesh fidelity analysis and refinement control
- **getStats.py**: Mesh quality metrics (aspect ratio, orthogonal quality, etc.)

### **postUtilities/** — Results Processing

- **pvPost.py**: Automated ParaView field extraction (pressure, velocity, etc.)
- **plotForces.py**: Force coefficient time-history and convergence plots
- **summary.py**: Compact results summary with aerodynamic coefficients
- **pptGeneration.py**: Automated PowerPoint report generation
- **postProReportGen.py**: LaTeX-based technical report generation

### **rideHeightUtils/** — Suspension Kinematics **[CORE NEW FEATURE]**

#### **doubleWishbonePushrodKinematics.py**
Kinematic solver for double-wishbone + pushrod suspension systems:
- Solves 3D linkage constraint equations
- Outputs per-ride-height: wheel center position, rotation (camber/toe), rocker angle, damper stroke
- Integrates with `utilities.py` for automatic geometry transformation

#### **suspensionHardpoints_template.cfg**
Template defining suspension geometry for each corner (FL, FR, RL, RR):

**Per-Corner Sections** `[SUSP_FL]`, `[SUSP_FR]`, `[SUSP_RL]`, `[SUSP_RR]`:
- **Kinematic Points**:
  - `WHEEL_CENTER_STATIC` — Static wheel position
  - `UCA_F_INNER`, `UCA_R_INNER`, `UCA_OUTER_STATIC` — Upper Control Arm chassis/wheel mounts
  - `LCA_F_INNER`, `LCA_R_INNER`, `LCA_OUTER_STATIC` — Lower Control Arm chassis/wheel mounts
  - `TIE_INNER`, `TIE_OUTER_STATIC` — Tie rod mounts
  - `PUSHROD_OUTER_STATIC` — Pushrod-to-rocker joint
  - `ROCKER_PIVOT`, `ROCKER_AXIS` — Rocker rotation center and axis
  - `ROCKER_PUSHROD_JOINT_REF`, `ROCKER_DAMPER_JOINT_REF`, `ROCKER_DAMPER_CHASSIS` — Rocker joint references
  - `WHEEL_AXIS_LOCAL`, `WHEEL_FORWARD_LOCAL` — Wheel orientation vectors

- **PID Keywords** (any key containing `PID`):
  - `FL_UCA_PID`, `FL_LCA_PID`, `FL_ROCKER_PID`, `FL_PUSHROD_PID`, `FL_DAMPER_PID`, `FL_TIE_PID`, `FL_WHEEL_PID`
  - Matched (case-insensitive, substring) against geometry file PIDs to apply component-specific transforms
  - Similar patterns for FR, RL, RR corners with fr/rl/rr prefixes
  - **Component classification** accepts the full name *and* common abbreviations, so values like `flprod` or `rlpr` are recognized as PUSHROD:

    | Component | Accepted keyword aliases |
    |-----------|--------------------------|
    | UCA       | `uca` |
    | LCA       | `lca` |
    | ROCKER    | `rocker`, `rock`, `rkr` |
    | PUSHROD   | `pushrod`, `prod`, `pshrd`, `push`, `pr` |
    | DAMPER    | `damper`, `damp`, `shock`, `strut` |
    | WHEEL     | `wheel`, `whl` |
    | TIE       | `tie`, `tierod`, `trod` |

---

## Setup & Configuration

### Case Directory Structure

```
CASES/
└── 001/
    ├── caseSetup                      # Configuration file (INI format)
    ├── constant/
    │   ├── triSurface/               # Geometry files (OBJ/STL)
    │   ├── transportProperties
    │   ├── turbulenceProperties
    │   └── ...
    ├── system/
    │   ├── blockMeshDict
    │   ├── snappyHexMeshDict
    │   ├── controlDict
    │   ├── fvSchemes
    │   ├── fvSolution
    │   └── ...
    ├── 0/                            # Initial conditions (created after mesh)
    └── postProcessing/               # (Generated after solver run)
```



### Ride-Height Input CSV

```csv
point,fl,fr,rl,rr,yaw,steer
1,0.05,0.05,0.02,0.02,0.0,0
2,0.10,0.10,0.05,0.05,0.0,5
3,-0.05,-0.05,-0.02,-0.02,0.0,0
```

Where:
- `fl, fr, rl, rr` = wheel vertical displacements (meters or mm depending on `RH_UNIT`)
- `yaw` = vehicle yaw angle (degrees, 0 for straight)
- `steer` = optional front road-wheel steering angle in degrees at the reference (front-left) wheel. The kinematic solver finds the rack travel that produces this angle at the front-left, then applies the same rack to the front-right so the Ackermann difference emerges. The tie-rod inner joint slides along the rack axis (`RACK_AXIS` in the hardpoint file, default lateral `0 1 0`) and the tie-rod length is preserved. Omit or set to `0` for no steer. Rear wheels are not steered.

---

## Per-Component Transform System **[NEW ARCHITECTURE]**

### How It Works

**Input**: Ride-height point (corner displacements + yaw angle)
↓
**Kinematic Solver**: Computes wheel position, rotation, rocker angle, damper stroke
↓
**Component Transforms**: For each suspension component, compute:
1. **Translation** — How the component moves
2. **Rotation** — Axis-angle rotation vector (radians)
3. **Pivot** — Center point for rotation
↓
**Geometry Transform**: Apply per-component transforms to matching PIDs in OBJ/STL files
↓
**Output**: Deformed geometry following suspension kinematics

### Component Transform Details

| Component | Translation | Rotation | Pivot |
|-----------|-------------|----------|-------|
| **WHEEL** | Wheel center displacement from solver | Rotation vector from solver (camber/toe) | Static wheel center |
| **UCA** | Outer point moves with wheel | Computed from linkage geometry change | Midpoint of chassis mounts |
| **LCA** | Outer point moves with wheel | Computed from linkage geometry change | Midpoint of chassis mounts |
| **ROCKER** | Zero (rotates in place) | Rocker angle × rocker axis | Rocker pivot |
| **PUSHROD** | Inherits from rocker | Inherits from rocker | Rocker pivot |
| **DAMPER** | Inherits from rocker | Inherits from rocker | Rocker pivot |
| **TIE** | Follows wheel translation | Computed from linkage geometry change | TIE_INNER (chassis mount) |

### Rotation Computation (UCA/LCA/TIE)

**Algorithm**: Rigid body fit using linkage vectors
1. **Static linkage vector**: From chassis mount to wheel attachment (static config)
2. **Current linkage vector**: Same, after wheel has moved
3. **Rotation axis**: Cross product of normalized vectors (perpendicular to both)
4. **Rotation angle**: Arc-cosine of dot product (angle between vectors)
5. **Rotation vector**: `axis × angle` (axis-angle representation in radians)

This ensures arms rotate realistically while maintaining distance constraints at their pivot points.

---

## Usage Workflow

### Step 1: Create New Case

```bash
cd CASES/001
python ../caseSetup.py --new
# Edit caseSetup file with your parameters
```

### Step 2: Enable Kinematic Solver (Optional)

Edit `caseSetup` and configure the `[RIDE_HEIGHT_SETUP]` block:
```ini
[RIDE_HEIGHT_SETUP]
RUN_RIDE_HEIGHT = True              # master switch for ride-height study
RIDE_HEIGHT_FILE = rideHeightConfig # CSV of ride-height points
RUN_RH_POINTS = 1 2 3 4 5           # subset of point indices to build
USE_DEPEND = True
RH_UNIT = mm                        # unit of corner displacements (mm | m)
INIT_RH =
RH_REF_WIDTH =
RH_REF_LEN =
USE_KINEMATIC_SOLVER = True         # enable double-wishbone + pushrod solver
HARDPOINT_FILE = rideHeightUtils/suspensionHardpoints_template.cfg
HARDPOINT_SCALE = 0.001             # 0.001 if hardpoints in mm, 1.0 if in m
SUSP_REQUIRE_ALL_PIDS_MATCHED = False
```

| Key | Description |
|-----|-------------|
| `RUN_RIDE_HEIGHT` | Master switch. `True` generates one child case per ride-height point |
| `RIDE_HEIGHT_FILE` | Name/path of the ride-height CSV file |
| `RUN_RH_POINTS` | Space-separated list of point indices (from the CSV) to build |
| `USE_DEPEND` | Reuse child-case dependencies/links when present |
| `RH_UNIT` | Unit of corner displacements in the CSV (`mm` or `m`) |
| `INIT_RH` | Optional initial/reference ride-height offset |
| `RH_REF_WIDTH`, `RH_REF_LEN` | Optional reference track width / wheelbase |
| `USE_KINEMATIC_SOLVER` | Enable per-component kinematic deformation (else simple wheel translation) |
| `HARDPOINT_FILE` | Path to the suspension hardpoint CFG (required if solver enabled) |
| `HARDPOINT_SCALE` | Coordinate-to-metre multiplier (`0.001` for mm, `1.0` for m) |
| `SUSP_REQUIRE_ALL_PIDS_MATCHED` | `True` errors on any unmatched part; `False` leaves it untouched with a warning |

Prepare hardpoint file (copy `suspensionHardpoints_template.cfg` and update coordinates):
```bash
cp rideHeightUtils/suspensionHardpoints_template.cfg constant/suspensionHardpoints.cfg
# Edit suspension geometry coordinates in CFG
```

Add PID keywords to geometry (any key containing `PID` in hardpoint CFG sections):
```ini
[SUSP_FL]
...
FL_WHEEL_PID = flwheel
FL_UCA_PID = fluca
FL_LCA_PID = fllca
FL_ROCKER_PID = flrocker
FL_PUSHROD_PID = flprod
FL_DAMPER_PID = fldamper
FL_TIE_PID = fltie
...
```
The keyword must be a **substring** of the actual part/group name in the OBJ/STL (e.g. keyword `fllca` matches a part named `body_susp_fllca`).

### Step 3: Setup Case

```bash
# Standard mesh + boundary condition setup. If RUN_RIDE_HEIGHT = True this
# automatically generates and builds the per-point child cases.
python ../caseSetup.py -s otr
```

> **Note:** `--rideHeightMode` is an internal flag that `caseSetup` passes to itself when re-running inside each child case. Do **not** invoke it manually — just set `RUN_RIDE_HEIGHT = True` and run `caseSetup` normally.

### Step 4: Process Ride-Height Points

If ride-height mode enabled, child cases created:
```
001_1/     # Case for ride-height point 1
001_2/     # Case for ride-height point 2
...
```

Each child case has:
- Deformed geometry following suspension kinematics
- Updated domain pitch/roll/heave
- Individual caseSetup with ride-height-specific parameters

### Step 5: Run Solver

```bash
cd 001_1/
blockMesh
snappyHexMesh
# Run solver (OpenFOAM command or cluster submission)
```

### Step 6: Post-Process Results

```bash
python ../postUtilities/postRun.py
# Generates:
# - Force convergence plots
# - Pressure/velocity field plots
# - Summary statistics
# - PowerPoint report (optional)
```

---

## Compilation

Use PyInstaller to create standalone executable (recommended for deployment):

```bash
bash compileFiles.sh
```

This creates `dist/caseSetup` executable with all templates and utilities bundled.

---

## Advanced Features

### Multi-Point Ride-Height Analysis

Generate cases for multiple ride-height points (e.g., acceleration, braking, cornering):

```csv
point,fl,fr,rl,rr,yaw
1,0.05,0.02,-0.02,-0.05,0.0      # Acceleration
2,-0.02,-0.05,0.05,0.02,0.0      # Braking
3,0.00,-0.08,0.08,0.00,5.0       # Right turn + roll
```

Run with `--rideHeightMode` (only for use by caseSetup itself internally!):
```bash
python caseSetup.py --rideHeightMode
```

Creates child cases for each point with automatically deformed suspension geometry.

### Mesh Fidelity Control

Define refinement regions in `defaultRefinements/`:
- Surface layers
- Local refinement zones
- Boundary layer y+ specification

### Custom Turbulence Models

Add new models in `zeroTemplates/models/`. Supported:
- k-omega SST
- Spalart-Allmaras (SA)
- k-epsilon (standard)

### Cluster Submission

Auto-generate SLURM/PBS scripts:
```bash
# Edit defaultCluster templates with your HPC settings
python caseSetup.py -s otr
# Generate: system/submitMesh.sh, system/submitSolver.sh, etc.
```

---

## Troubleshooting

### Kinematic Solver Fails
- Check hardpoint CFG syntax (INI format, exact key names)
- Verify all required keys present (see template)
- Ensure suspension geometry is physically feasible (no link crossing)

### Geometry Misaligned After Transform
- Verify PID keywords match geometry group names
- Check ride-height input (excessive corner displacements?)
- Review component pivot points in hardpoint CFG

### "PID scan skipped" / Unmatched PID
- "PID scan skipped" means a geometry could not be processed by the PID-aware path (unsupported or binary format); it falls back to simple wheel translation. Use **ASCII** OBJ/STL for suspension parts.
- An "Unmatched PID" warning means a hardpoint keyword does not appear in any part/group name. Confirm the keyword is a substring of the actual part name, or add an alias.
- Inspect the per-part transform log written next to each output geometry (`*.susp_pid_transform.json` / `.csv`) to confirm which parts were moved.
- Verify `HARDPOINT_SCALE` matches the unit of your hardpoint coordinates (`0.001` for mm).

### Memory Issues with Large Meshes
- Enable parallel decomposition in `system/decomposeParDict`
- Reduce mesh refinement in `snappyHexMeshDict`
- Split case into multiple smaller domains

---

## Performance & Scalability

- **Mesh generation**: ~minutes (depends on refinement)
- **Kinematic solve**: ~milliseconds per ride-height point
- **Geometry transform**: ~seconds per component (OBJ/STL file size dependent)
- **CFD solver**: Hours to days (HPC recommended for production)

---

## License & Attribution

This framework builds on OpenFOAM and community contributions. See individual module headers for specific acknowledgments.

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| **4.2-dev** | 2026 | Per-component transforms, wheel camber/toe, UCA/LCA/TIE rotation computation; PID keyword abbreviation matching (e.g. `flprod`/`rlpr` → PUSHROD); macOS multiprocessing fix (`__main__` guard) preventing child-case process recursion; PID transform log `pivot=None` fix |
| 4.1.1 | 2025 | Kinematic solver integration, hardpoint CFG support |
| 4.0 | 2024 | Initial multi-setup template framework |

---

## Support & Contact

For issues, questions, or contributions, refer to project documentation or contact maintainers.

