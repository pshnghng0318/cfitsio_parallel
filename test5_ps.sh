#!/bin/bash

USE_SUBSET=0  # 0: fits_reader_mpi; 1: fits_reader_mpi_subset
FITSFILE=/images/set_QA/alma_band1_orionkl_bandscan_combined.fits
X=10  # Number of processes

if [[ $USE_SUBSET -eq 1 ]]; then
    EXEC="fits_reader_mpi_subset"
    SRC="fits_reader_mpi_subset.cpp"
    RESULT_FILE="result${X}_subset.dat"
    PLOT="plotresults_mpi_subset.py"
else
    EXEC="fits_reader_mpi"
    SRC="fits_reader_mpi.cpp"
    RESULT_FILE="result${X}.dat"
    PLOT="plotresults_mpi.py"
fi

# Compile
/opt/homebrew/bin/g++-15 -std=c++20 \
    -I/opt/homebrew/Cellar/open-mpi/5.0.7/include \
    -lmpi -L/opt/homebrew/Cellar/open-mpi/5.0.7/lib \
    -O3 -ltbb -L/opt/homebrew/Cellar/tbb/2022.1.0/lib/ -I/opt/homebrew/Cellar/tbb/2022.1.0/include \
    -o $EXEC $SRC -lcfitsio \
    -L/opt/homebrew/Cellar/cfitsio/4.6.2/lib \
    -I/opt/homebrew/Cellar/cfitsio/4.6.2/include -pthread

# If compilation successful
if [ $? -eq 0 ]; then
    > "$RESULT_FILE"  # Clear previous results

    for i in {1..5}; do
        echo "========== Iteration $i =========="

        sudo purge
        sync && sleep 1

        MEM_LOG="mem_iter${i}_${X}.log"
        MEM_PNG="mem_iter${i}_${X}.png"

        # Use psrecord to monitor the memory usage
        psrecord "mpirun -np $X ./$EXEC $FITSFILE" \
            --interval 1 \
            --include-children \
            --log "$MEM_LOG" \
            --plot "$MEM_PNG"

        echo "Saved memory usage to $MEM_LOG and $MEM_PNG"
    done

    # Run plotting scripts (in background)
    python plotspectrum.py &
    python "$PLOT" "$X" &
    python plotratio.py "$X" &
fi
