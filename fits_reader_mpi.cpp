// mpic++ -std=c++11 -o fits_reader_mpi fits_reader_mpi.cpp -lcfitsio -L/opt/homebrew/Cellar/cfitsio/4.6.2/lib -I/opt/homebrew/Cellar/cfitsio/4.6.2/include -pthread
#include <mpi.h>
#include <iostream>
#include <vector>
#include <fitsio.h>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <execution>
#include <cmath>
#include <sys/types.h>
#include <sys/sysctl.h>

void read_fits_chunk(const char* filename, long start_pixel, long num_pixels, std::vector<float>& buffer) {

    int status = 0;
    //float nulval = 0.0;

    fitsfile* fptr;

    if (fits_open_file(&fptr, filename, READONLY, &status)) {
        fits_report_error(stderr, status);
        MPI_Abort(MPI_COMM_WORLD, status);
    }
    // +1: fits start pixel is 1-based
    if (fits_read_img(fptr, TFLOAT, start_pixel + 1, num_pixels, nullptr, buffer.data(), nullptr, &status)) {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        MPI_Abort(MPI_COMM_WORLD, status);
    }
    //fits_read_subset(fptr, TFLOAT, start_pixel + 1, num_pixels, &nulval, buffer.data(), NULL, &status);
    //fits_read_pix(fptr, TFLOAT, start_pixel + 1, num_pixels, &nulval, buffer.data(), NULL, &status);
    
    fits_close_file(fptr, &status);
}

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    long i;

    // Get available CPU number
    int rank = 0, size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Use half memory size
    int64_t mem_size = 0;
    size_t len = sizeof(mem_size);
    sysctlbyname("hw.memsize", &mem_size, &len, nullptr, 0);

    if (rank == 0) {
        std::cout << "Total: " << size << " Perf CPU and " << mem_size / (1024*1024*1024) << " GB Memory" << std::endl;
    }
    mem_size = mem_size * 1 / 8 ;
    //mem_size = mem_size * 1 / 4 ;
    //mem_size /= 2;
    //mem_size = mem_size * 3 / 4 ;

    if (rank == 0) {
        std::cout << "Using half memory size: " << mem_size / (1024*1024*1024) << "GB" << std::endl;
    }

    std::ofstream outfile("result" + std::to_string(size) + ".dat", std::ios::app);

    if (argc != 2) {
        if (rank == 0) {
            std::cerr << "Usage: mpirun -np <num_procs> ./fits_reader_mpi <FITS file>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    const char* filename = argv[1];
    fitsfile* fptr;
    int status = 0, bitpix = 0, naxis = 0;
    long naxes[4] = {1, 1, 1, 1};

    std::int64_t npixels = 1;
    int n_channels = 1, n_x = 1, n_y = 1;
    std::int64_t less_pixels = 1;
    std::int64_t loop_pixels = 1;
    int less_channels = 1;
    int loop_channels = 1;
    //std::cout << "npixels: " << npixels << ", n_x: " << n_x << ", n_y: " << n_y << ", n_channels: " << n_channels << std::endl;

    if (rank == 0) {
        if (fits_open_file(&fptr, filename, READONLY, &status)) {
            fits_report_error(stderr, status);
            MPI_Abort(MPI_COMM_WORLD, status);
        }
        if (fits_get_img_param (fptr, 4, &bitpix, &naxis, naxes, &status)) {
            fits_report_error(stderr, status);
            fits_close_file(fptr, &status);
            MPI_Abort(MPI_COMM_WORLD, status);
        }
        fits_close_file(fptr, &status);

        n_x = naxes[0];
        n_y = naxes[1];
        n_channels = (naxis >= 3) ? naxes[2] : 1;
        npixels = std::int64_t(n_x) * n_y * n_channels;

        int nbytes = 4;
        switch (bitpix) {
            case 8: nbytes = 1; break;   // BYTE
            case 16: nbytes = 2; break;  // SHORT
            case 32: nbytes = 4; break;  // LONG
            case -32: nbytes = 4; break; // FLOAT
            case -64: nbytes = 8; break; // DOUBLE
            default:
                if (rank == 0) {
                    std::cerr << "Unsupported BITPIX value: " << bitpix << std::endl;
                }
                MPI_Finalize();
                return -1; 
        }
        
        if (npixels * nbytes > mem_size) {
            less_pixels = mem_size / nbytes;
            less_channels = less_pixels / (n_x * n_y);
        } else {
            less_pixels = npixels;
            less_channels = n_channels;
        }

        std::cout << "Total pixels: " << npixels << ", n_x: " << n_x << ", n_y: " << n_y 
              << ", n_channels: " << n_channels << ", nbytes: " << nbytes 
              << ", less_pixels: " << less_pixels << std::endl;

    }

    // Broadcast metadata
    MPI_Bcast(&npixels, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_y, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_channels, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bitpix, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&less_pixels, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&less_channels, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    double io_time = 0.0;
    double compute_time = 0.0;
    double total_time = 0.0;
    MPI_Bcast(&io_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&compute_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&total_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    long total_loops = std::int64_t(npixels) / less_pixels + ((std::int64_t(npixels) % less_pixels) ? 1 : 0);
    if (rank == 0) {
        std::cout << "Total loops: " << total_loops << std::endl;
    }
    std::vector<float> my_buffer;
    std::vector<double> local_spectrum(n_channels, 0.0);
    std::vector<long> valid_pixels(n_channels, 0);


    for (long nloop = 0; nloop < total_loops; ++nloop) {
        if (nloop == total_loops - 1) {
            loop_pixels = npixels - nloop * less_pixels;
        } else {
            loop_pixels = less_pixels;
        }
        loop_channels = loop_pixels / (n_x * n_y);
        long base_chunk = loop_pixels / size; // pixels per rank
        long remainder = loop_pixels % size;
        long channels_per_rank = loop_channels / size;
        long extra = loop_channels % size;
        if (rank == 0) {
            std::cout << "base_chunk: " << base_chunk 
                      << ", remainder: " << remainder 
                      << ", channels_per_rank: " << channels_per_rank 
                      << ", extra: " << extra 
                      << ", loop_channels: " << loop_channels 
                      << std::endl;
        }

        long my_first_channel = rank * channels_per_rank + std::min((long)rank, extra) + nloop * loop_channels;
        long my_last_channel = my_first_channel + channels_per_rank - 1;
        //std::cout << "my_first_channel: " << my_first_channel << ", my_last_channel: " << my_last_channel << ", rank: " << rank << ", size: " << size << std::endl;
        if (rank < extra) {
            my_last_channel += 1;
        }
        if (my_last_channel >= n_channels) {
            my_last_channel = n_channels - 1; // Ensure we don't exceed the total number of channels
        }

        long my_channels = my_last_channel - my_first_channel + 1;
        long my_pixels = my_channels * n_x * n_y;
        long my_start_pixel = my_first_channel * n_x * n_y;
        my_buffer.resize(my_pixels, 0.0f);
        my_buffer.shrink_to_fit();


        //std::cout << "rank: " << rank << ", my_channels: " << my_channels << ", my_pixels: " << my_pixels << ", my_start: " << my_start_pixel << std::endl;

        // Start timing (total time)
        MPI_Barrier(MPI_COMM_WORLD);
        auto total_start = std::chrono::high_resolution_clock::now();

        // Time I/O
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processes before starting I/O
        auto io_start = std::chrono::high_resolution_clock::now();

        read_fits_chunk(filename, my_start_pixel, my_pixels, my_buffer);

        MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processes after I/O
        auto io_end = std::chrono::high_resolution_clock::now();
        if (rank == 0) {
            io_time += std::chrono::duration<double>(io_end - io_start).count();
            std::cout << io_time << std::endl;
        }
    
        // Time computation
        MPI_Barrier(MPI_COMM_WORLD);        
        auto compute_start = std::chrono::high_resolution_clock::now();

        // local spectrum

        using std::isnan;
        for (long iter = 0; iter < my_pixels; ++iter) {
            if (!isnan(my_buffer[iter])) {
                long channel_index = (my_start_pixel + iter) / (n_x * n_y);
                local_spectrum[channel_index] += my_buffer[iter];
                valid_pixels[channel_index] += 1;
            }
        }
        my_buffer.clear();

        MPI_Barrier(MPI_COMM_WORLD);
        auto compute_end = std::chrono::high_resolution_clock::now();
        auto total_end = std::chrono::high_resolution_clock::now();

        if (rank == 0) {
            compute_time += std::chrono::duration<double>(compute_end - compute_start).count();
            total_time += std::chrono::duration<double>(total_end - total_start).count();
        }
    }

    // Time computation
    MPI_Barrier(MPI_COMM_WORLD);        
    auto compute_start = std::chrono::high_resolution_clock::now();
    auto total_start = std::chrono::high_resolution_clock::now();
    //// Normalize the local spectrum
    //for (long c = 0; c < n_channels; ++c) {
    //    if (valid_pixels[c] > 0) {
    //        local_spectrum[c] /= valid_pixels[c];
    //    }
    //}

    // Normalize the local spectrum using std::transform
    std::transform(std::execution::par, 
        local_spectrum.begin(), 
        local_spectrum.end(), 
        valid_pixels.begin(), 
        local_spectrum.begin(), 
        [](float spectrum, int count) {
            return (count > 0) ? spectrum / count : spectrum;
        }
    );

    MPI_Barrier(MPI_COMM_WORLD);

    // Reduce the local spectrum to a global spectrum
    std::vector<double> global_spectrum(n_channels, 0.0); // Initialize a global spectrum
    // Assign rank 0 to reduce the local spectrum
    MPI_Reduce(local_spectrum.data(), global_spectrum.data(), n_channels, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 

    MPI_Barrier(MPI_COMM_WORLD);
    auto compute_end = std::chrono::high_resolution_clock::now();
    auto total_end = std::chrono::high_resolution_clock::now();

    if (rank == 0) {
        compute_time += std::chrono::duration<double>(compute_end - compute_start).count();
        total_time += std::chrono::duration<double>(total_end - total_start).count();
    }
    if (rank == 0) {
        std::cout << "I/O time:      " << io_time << " seconds\n";
        std::cout << "Compute time:  " << compute_time << " seconds\n";
        std::cout << "Total runtime: " << total_time << " seconds\n";
        outfile << total_time << " " << io_time << " " << compute_time << " sec." << std::endl;
    }
    outfile.close();

    if (rank == 0) {
        std::ofstream specout("spectrum_mpi.dat");
        for (long c = 0; c < n_channels; ++c) {
            specout << c << "\t" << global_spectrum[c] << "\n";
        }
        specout.close();
        //std::cout << "Spectrum written to spectrum_mpi.dat\n";
    }

    MPI_Finalize();

    return 0;
}
