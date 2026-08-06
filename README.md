# Adaptive Pn-wave Tomography Workflow

A workflow for uppermost mantle Pn-wave velocity and anisotropy tomography using adaptive non-uniform grids.

## Overview

This repository provides preprocessing, postprocessing and visualization tools developed for Pn-wave travel-time tomography.

The workflow implements an adaptive two-level non-uniform grid strategy, where grid resolution is determined according to seismic ray density. This approach improves model resolution in regions with dense ray coverage while maintaining inversion stability in poorly sampled areas.

Example result in Eastern Tibet:
![alt text](<Figure 8.jpg>)

## Features

- Adaptive two-level non-uniform grid construction
- Pn phase data preprocessing and quality control
- Epicentral distance calculation and ray selection
- Pn velocity and anisotropy model visualization
- Checkerboard resolution test generation
- GMT-based automatic plotting workflow

## Workflow

Input phase data

↓

Phase preprocessing (`rdphase`)

↓

Pn tomography inversion (`ctomo`, not included)

↓

Solution extraction (`readsol`)

↓

GMT visualization

## Repository Structure

src/
rdphase.c
readsol.c
syphase.c

scripts/
rdphase.sh
ctomo.sh
read.sh
sol.sh
ani.sh
sks.sh

figures/

examples/

## Requirements

- GCC compiler
- GMT 6.x
- Bash environment
- Linux system recommended

## Usage

### 1. Phase preprocessing

Modify parameters in: rdphase.sh

including:

- geographic range
- input phase file
- filtering criteria

Run:

```bash
./rdphase.sh
```

### 2. Tomography inversion

The inversion module ctomo performs:

adaptive grid generation
ray tracing
damped inversion

The source code is not included in this repository.

Run:

```bash
./ctomo.sh
```

### 3. Result processing

Run:

```bash
./read.sh
```

to generate:

sol.dat &
ani.txt

### 4. Visualization

Generate velocity and anisotropy maps:

```bash
./sol.sh
./ani.sh
```

### Resolution Test

```bash
syphase.c
```

 generates checkerboard models based on the actual earthquake-station ray distribution for resolution assessment.

### Applications

The workflow has been applied to investigate uppermost mantle structures beneath continental collision zones, including:

Eastern Mediterranean-Anatolia region:

Uppermost mantle nonuniform grid Pn velocity and anisotropy tomography beneath the Eastern Mediterranean and Anatolian region[J]. Earth and Planetary Science Letters, 2026, 690, 120166. DOI: 10.1016/j.epsl.2026.120166

Tibetan Plateau:

Pn velocity and anisotropy tomography with nonuniform grid beneath the Tibetan Plateau and adjacent regions[J]. Tectonophysics, 2026, 935, 231278. DOI: 10.1016/j.tecto.2026.231278
