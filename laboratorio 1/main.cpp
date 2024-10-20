#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cmath> 
#include <omp.h>
#include <fitsio.h>
#include <vector>

using namespace std;

// Data functions
float ox(int tx, int ux) { return (tx + ux) / 2.0f; }
float ux(int ty, int uy) { return (ty + uy) / 2.0f; }
float alpha(int tx, int ty, int ux, int uy) { return sqrt(pow((ux - tx), 2) + pow((uy - ty), 2)) / 2.0f; }
float theta(int tx, int ty, int ux, int uy) { return atan2((uy - ty), (ux - tx)); }

int main(int argc, char **argv) {

    // terminal input
    int opt;
    char* filename = nullptr;
    int alpha = 0, percentage = 0, betha = 0, levelOneThreads = 0, levelTwoThreads = 0;

    // get options
    while ((opt = getopt(argc, argv, "f:a:r:b:u:d:")) != -1) {
        switch (opt) {
            case 'f':
                filename = optarg;
                break;
            case 'a':
                alpha = std::atoi(optarg);
                break;
            case 'r':
                percentage = std::atoi(optarg);
                break;
            case 'b':
                betha = std::atoi(optarg);
                break;
            case 'u':
                levelOneThreads = std::atoi(optarg);
                break;
            case 'd':
                levelTwoThreads = std::atoi(optarg);
                break;
            case '?':
                std::cerr << "Unknown option: " << char(optopt) << std::endl;
                return 1;
            default:
                std::cerr << "Use: " << argv[0] << " -f <filename> -a <minimun alpha> -r <relative vote> -b <number of bethas> -u <number of threads 1> -d <number of threads 2>" << std::endl;
        }
    }

    if (filename == nullptr) {
        cerr << "Error: Debe proporcionar un nombre de archivo con -f." << endl;
        return 1;
    }

    // image variables
    fitsfile *fptr;
    int status = 0;
    long fpixel = 1, naxis = 2;
    long naxes[2];

    // Open fits file
    if (fits_open_file(&fptr, filename, READONLY, &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    // Get image size
    if (fits_get_img_size(fptr, 2, naxes, &status)) {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return status;
    }

    // Memory for image
    double *image = (double *) malloc(naxes[0] * naxes[1] * sizeof(double));
    if (image == nullptr) {
        cerr << "Error al asignar memoria para la imagen." << endl;
        fits_close_file(fptr, &status);
        return 1;
    }
 
    // Vector that will contain every index of a 
    vector<long> pixels;

    #pragma omp parallel num_threads(levelOneThreads) shared(image, naxes, fptr, pixels)
    {
        int thread_id = omp_get_thread_num();
        long totalPixels = naxes[0] * naxes[1];
        long pixelsPerThread = (totalPixels + levelOneThreads - 1) / levelOneThreads;
        long start = thread_id * pixelsPerThread + 1;
        long end = min(start + pixelsPerThread - 1, totalPixels);
        long pixelsRead = end - start + 1;

        int local_status = 0; 
        fits_read_img(fptr, TDOUBLE, start, pixelsRead, NULL, &image[start - 1], NULL, &local_status);

        vector<long> localPixels;

        for(int i = start - 1; i < end; i++){
            if(image[i] == 255){
                localPixels.push_back(i);
            }
        }

        #pragma omp critical
        pixels.insert(pixels.end(), localPixels.begin(), localPixels.end());

    }
    
    cout << pixels.size() << endl;

    // Cerrar el archivo FITS
    if (fits_close_file(fptr, &status)) {
        fits_report_error(stderr, status);
        free(image);
        return status;
    }

    // Imprimir cada pixel de la imagen
    double num = 0;
    #pragma omp parallel for num_threads(levelOneThreads) shared(num)
    for (long i = 0; i < naxes[0] * naxes[1]; i++) {
        #pragma omp critical
        if(image[i] == 255){
            // cout << "Pixel " << i << ": " << image[i] << endl;    
            num = num + 1;
        }
    }
    cout << num << endl;

    // Liberar la memoria
    free(image);

    return 0;
}
