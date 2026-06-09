#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e perpendicular.hdf5 ]
then
    python3 makeIC.py
fi


../../../swift --hydro --limiter --sync --self-gravity --black-holes --threads=16 --pin  perpendicular.yml

