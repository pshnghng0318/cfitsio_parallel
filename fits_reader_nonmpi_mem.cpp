// g++ -std=c++11 -o fits_reader_nonmpi fits_reader_nonmpi.cpp -I/opt/cfitsio_non_reentrant/include -L/opt/cfitsio_non_reentrant/lib
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
    }
    // +1: fits start pixel is 1-based 
    if (fits_read_img(fptr, TFLOAT, start_pixel + 1, num_pixels, nullptr, buffer.data(), nullptr, &status)) {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
    }

    // Replace NaN with 0.0
    //std::replace_if(buffer.begin(), buffer.end(), [](float val){ return std::isnan(val); }, 0.0f); 
    
    fits_close_file(fptr, &status);

}

int main(int argc, char* argv[]) {
    long i;

    int64_t mem_size = 0;
    size_t len = sizeof(mem_size);
    sysctlbyname("hw.memsize", &mem_size, &len, nullptr, 0);
    
    std::ofstream outfile("result1.dat", std::ios::app);

    if (argc != 2) {
        std::cerr << "Usage: g++ ./fits_reader_nonmpi <FITS file>" << std::endl;
        return 1;
    }
    mem_size = mem_size * 1 / 4 ;

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

    if (fits_open_file(&fptr, filename, READONLY, &status)) {
        fits_report_error(stderr, status);
    }
    if (fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status)) {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
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
            std::cerr << "Unsupported BITPIX value: " << bitpix << std::endl;
            return -1; 
    }
        
    if (npixels * nbytes > mem_size) {
        less_pixels = mem_size / nbytes;
        less_channels = less_pixels / (n_x * n_y);
    } else {
        less_pixels = npixels;
        less_channels = n_channels;
    }
    double io_time = 0.0;
    double compute_time = 0.0;
    double total_time = 0.0;
    long total_loops = std::int64_t(npixels) / less_pixels + ((std::int64_t(npixels) % less_pixels) ? 1 : 0);
    std::vector<float> buffer;
    std::vector<double> spectrum(n_channels, 0.0);
    std::vector<long> valid_pixels(n_channels, 0);

    for (long nloop = 0; nloop < total_loops; ++nloop) {
        if (nloop == total_loops - 1) {
            loop_pixels = npixels - nloop * less_pixels;
        } else {
            loop_pixels = less_pixels;
        }
        loop_channels = loop_pixels / (n_x * n_y);
        long base_chunk = loop_pixels; // pixels per rank
        long channels_per_rank = loop_channels;
        long my_first_channel = nloop * loop_channels;
        long my_last_channel = my_first_channel + channels_per_rank - 1;
        long my_first_pixel = my_first_channel * n_x * n_y;
        long my_last_pixel = my_last_channel * n_x * n_y + n_x * n_y - 1;

        long my_channels = my_last_channel - my_first_channel + 1;
        long my_pixels = my_channels * n_x * n_y;
        long my_start_pixel = my_first_channel * n_x * n_y;
        buffer.resize(my_pixels, 0.0f);
        buffer.shrink_to_fit();

        std::vector<float> buffer(npixels);
        // Start timing (total time)
        auto total_start = std::chrono::high_resolution_clock::now();
        
        // Time I/O
        auto io_start = std::chrono::high_resolution_clock::now();
        read_fits_chunk(filename, my_start_pixel, my_pixels, buffer);

        auto io_end = std::chrono::high_resolution_clock::now();
        io_time += std::chrono::duration<double>(io_end - io_start).count();

        // Time computation
        auto compute_start = std::chrono::high_resolution_clock::now();

        
        using std::isnan;
        for (long iter = 0; iter < my_pixels; ++iter) {
            if (!isnan(buffer[iter])) {
                long ch = (my_start_pixel+iter) / (n_x * n_y);
                spectrum[ch] += buffer[iter];
                valid_pixels[ch] += 1;
            }
        }
        buffer.clear();

        auto compute_end = std::chrono::high_resolution_clock::now();
        auto total_end = std::chrono::high_resolution_clock::now();

        compute_time += std::chrono::duration<double>(compute_end - compute_start).count();
        total_time += std::chrono::duration<double>(total_end - total_start).count();
    }

    auto compute_start = std::chrono::high_resolution_clock::now();
    auto total_start = std::chrono::high_resolution_clock::now();

    std::transform(std::execution::par, 
        spectrum.begin(), 
        spectrum.end(), 
        valid_pixels.begin(), 
        spectrum.begin(), 
        [](float spectrum, int count) {
            return (count > 0) ? spectrum / count : spectrum;
        }
    );

    auto compute_end = std::chrono::high_resolution_clock::now();
    auto total_end = std::chrono::high_resolution_clock::now();

    // Report timings
    compute_time += std::chrono::duration<double>(compute_end - compute_start).count();
    total_time += std::chrono::duration<double>(total_end - total_start).count();

    std::cout << "I/O time:      " << io_time << " seconds\n";
    std::cout << "Compute time:  " << compute_time << " seconds\n";
    std::cout << "Total runtime: " << total_time << " seconds\n";
    outfile << total_time << " " << io_time << " " << compute_time << " sec." << std::endl;
    outfile.close();

    std::ofstream specout("spectrum_nonmpi.dat");
    for (long c = 0; c < n_channels; ++c) {
        specout << c << "\t" << spectrum[c] << "\n";
    }
    specout.close();
    //std::cout << "Spectrum written to spectrum_nonmpi.dat\n";    
    
    return 0;
}
