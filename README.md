# Adaptive Pn-wave Tomography Workflow

A workflow for uppermost mantle Pn-wave velocity and anisotropy tomography using adaptive non-uniform grids.

## Overview

This repository provides preprocessing, postprocessing and visualization tools developed for Pn-wave travel-time tomography.

The workflow implements an adaptive two-level non-uniform grid strategy, where grid resolution is determined according to seismic ray density. This approach improves model resolution in regions with dense ray coverage while maintaining inversion stability in poorly sampled areas.

Example result in Eastern Tibet:
![alt text](<figures/Result in Eastern Tibet.jpg>)

## Features

- Adaptive two-level non-uniform grid construction
- Pn phase data preprocessing and quality control
- Epicentral distance calculation and ray selection
- Pn velocity and anisotropy model visualization
- Checkerboard resolution test generation
- GMT-based automatic plotting workflow

Schematic diagram of ray tracing with two layers of non-uniform grid:
![alt text](<figures/Ray tracing.png>)
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

**src/**
rdphase.c
readsol.c
syphase.c

**scripts/**
rdphase.sh
ctomo.sh
read.sh
sol.sh
ani.sh
sks.sh

**figures/**

**examples/**

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

**The source code is not included in this repository.** You can obtain the core inversion code from **the SI of the two articles** mentioned in Applications.

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
 
 Different synthesis models and related parameter settings are explained in the program comments.

## C-CODES
### rdphase.c:
	Sorts the data and outputs binary data structures stations.b, events.b & phases.b.  There are a bunch of parameters at the start of this code which define the raypath (usually Pn).  Make sure that only refracted raypaths of the same phase are used.  The most important parameters are STNMIN & EVNMIN which set the minimum arrivals required at each station and event.  I prefer these be at least 10.  Ath the end of rdphase.c is a list of event selection criteria that the user should check out.  These often get modified.  The ISC data, for example, had some quality criteria that I used to select data.  Other data sets may have differing criteria.
	At present rdphase is configured to read data files in the *.b format as defined in ./src/data14.h.  Note that stations.b and events.b have the same format, it's just that stations.b has info about each station and events.b have info about each event.
	Input used stdin, output is stations.b events.b & phases.b.
	Only flag is "-region xmin/xmax/ymin/ymax" where xmin, xmax, ymin, & ymax are the longitude latitude bounds.

### syphase.c:
	This program creates synthetic data structures.  It is not well written.  It creates checkerboard patterns at present.  You will need to hack at it to change the test pattern.
	Flags are:
		"-region xmin/xmax/ymin/ymax" where xmin, xmax, ymin, & ymax are the longitude latitude bounds.
		"-mesh xmesh/ymesh" where xmesh & ymesh define the grid size in cells/degree.
		"-cellmesh cmesh" where cell mesh is the test square size in squares/degree.
		
### readpha.c:
	Reads phases.b and outputs.  phases.b contains info for each arrival in the data set.  It uses stdout.  Use the output to plot raypaths.
	Flags are:
		"-file filename" filename is input filename.  A default is hard coded in.
		"-rays" sets output format for input to psxy to draw raypaths.
	
### readsta.c:  
	Reads stations.b or events.b from the file specified by the -file flag.  Use this to output for plotting station delays.
	Flags are:
		"-file filename" filename is input filename.  A default is hard coded in.

### readsol.c:  
	Reads sol.b from the file specified by the -file flag.  Use this to plot the slowness image.
	Flags are:
		"-file filename" filename is input filename.  A default is hard coded in.
		"-region xmin/xmax/ymin/ymax" where xmin, xmax, ymin, & ymax are the longitude latitude bounds.
		"-mesh xmesh/ymesh" where xmesh & ymesh define the grid size in cells/degree.
		"-vn velocity" where velocity is the Pn velocity (to convert slowness perturbations to velocity purturbations.


## Applications

The workflow has been applied to investigate uppermost mantle structures beneath continental collision zones, including:

**Eastern Mediterranean-Anatolia region**:

Uppermost mantle nonuniform grid Pn velocity and anisotropy tomography beneath the Eastern Mediterranean and Anatolian region[J]. Earth and Planetary Science Letters, 2026, 690, 120166. DOI: 10.1016/j.epsl.2026.120166
![alt text](<figures/Mediterranean graphical abstract.png>)

**Tibetan Plateau**:

Pn velocity and anisotropy tomography with nonuniform grid beneath the Tibetan Plateau and adjacent regions[J]. Tectonophysics, 2026, 935, 231278. DOI: 10.1016/j.tecto.2026.231278
![alt text](<figures/Tibet graphical abstract.png>)