#!/bin/bash

FITSFILE=/images/set_QA/alma_band1_orionkl_bandscan_combined.fits

# clear for appending
> "result1.dat"

/opt/homebrew/bin/g++-15 -std=c++20 \
    -O3 -ltbb -L/opt/homebrew/Cellar/tbb/2022.1.0/lib/ -I/opt/homebrew/Cellar/tbb/2022.1.0/include \
    -o fits_reader_nonmpi_mem fits_reader_nonmpi_mem.cpp -I/opt/cfitsio_non_reentrant/include -L/opt/cfitsio_non_reentrant/lib -lcfitsio

#g++ -std=c++17 -o fits_reader_nonmpi_mem fits_reader_nonmpi_mem.cpp -I/opt/cfitsio_non_reentrant/include -L/opt/cfitsio_non_reentrant/lib -lcfitsio

# non-reentrant cfitsio
for i in {1..5}; do
    sudo purge
    sync && sleep 1

    MEM_LOG="mem_1_iter_${i}.log"
    MEM_PNG="mem_1_iter_${i}.png"

    # Use psrecord to monitor the memory usage
    psrecord "DYLD_LIBRARY_PATH=/opt/cfitsio_non_reentrant/lib ./fits_reader_nonmpi_mem $FITSFILE" \
        --interval 1 \
        --include-children \
        --log "$MEM_LOG" \
        --plot "$MEM_PNG"
done

#python plotspectrum.py &
python plotresults_mpi.py 1 &
#python plotratio.py 1
