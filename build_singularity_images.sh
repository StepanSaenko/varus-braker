#!/bin/bash

singularity build braker3.sif docker://katharinahoff/bioinformatics-notebook:devel # contains BRAKER3 + OMArk and other software
singularity build tetools.sif docker://dfam/tetools:latest # RepeatModeler2 + RepeatMasker
singularity build busco.sif docker://quay.io/biocontainers/busco:latest # BUSCO
# singularity build varus.sif docker://katharinahoff/varus-notebook:devel # not functional, yet
