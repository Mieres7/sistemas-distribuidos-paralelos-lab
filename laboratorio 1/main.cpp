#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cmath> 
#include <omp.h>
#include <fitsio.h>
#include <vector>

using namespace std;


typedef struct Ellipse {
    float ox;
    float xy;
    float alpha;
    float betha;
    float tetha;
};

// Data functions
float getOx(int tx, int ux) { return (tx + ux) / 2.0f; }
float getOy(int ty, int uy) { return (ty + uy) / 2.0f; }
float getAlpha(int tx, int ty, int ux, int uy) { return (sqrt(pow((ux - tx), 2) + pow((uy - ty), 2)) / 2.0f); }
float getTheta(int tx, int ty, int ux, int uy) { return atan2((uy - ty),(ux - tx)); }

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

    fits_read_img(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1], NULL, image, NULL, &status);
 
    // Vector that will contain every pixel with a value o 255, x and y are coordinates
    vector<pair<int,int>> pixels;

    // Vector with found ellipses
    vector<Ellipse> ellipses;

    // First parallel level - create edge list 
    omp_set_nested(0);
    #pragma omp parallel num_threads(levelOneThreads) shared(image, naxes, fptr, pixels)
    {   
        // extraction of pixels with value 255, generating its coordinates
        vector<pair<int,int>> localPixels;
        #pragma omp for
        for(int i = 0; i < naxes[0]*naxes[1]; i++){
            if(image[i] == 255){
                int x = i % naxes[0];
                int y = i / naxes[0];
                pair<int, int> pixel;
                pixel.first = x;
                pixel.second = y;
                localPixels.push_back(pixel);
            } 
        }
        #pragma omp critical // insert every pixel in vector pixels
        pixels.insert(pixels.end(), localPixels.begin(), localPixels.end());

        // Second parallel level - Hough algorithm
        #pragma omp parallel num_threads(levelTwoThreads) shared(pixels, ellipses)
        {
            #pragma omp for
            for(pair<int,int> pixel_t: pixels ){

                #pragma omp for
                for(pair<int, int> pixel_u : pixels){

                    int vote[pixels.size()]; // array for betha votes

                    float oX = getOx(pixel_t.first, pixel_u.first);
                    float xY = getOy(pixel_t.second, pixel_u.second);
                    float alpha = getAlpha(pixel_t.first, pixel_t.second, pixel_u.first, pixel_u.second);
                    float theta = getTheta(pixel_t.first, pixel_t.second, pixel_u.first, pixel_u.second);

                    #pragma omp for
                    for(pair<int, int> pixel_k : pixels){
                        
                        if(pixel_k == pixel_t && pixel_k == pixel_u){
                            continue;
                        }

                        


                    }


                }

            }
        }




    }
    
    cout << pixels.size() << endl;

    // Close fits file
    if (fits_close_file(fptr, &status)) {
        fits_report_error(stderr, status);
        free(image);
        return status;
    }

    free(image);

    return 0;
}
