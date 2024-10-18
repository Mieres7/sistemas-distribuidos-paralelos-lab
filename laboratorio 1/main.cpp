#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cmath> 

#include <fitsio.h>

using namespace std;


// Data functions
float ox(int tx, int ux){return (tx + ux)/2;}
float ux(int ty, int uy){return (ty + uy)/2;}
float alpha(int tx, int ty, int ux, int uy){return sqrt(pow((ux - tx), 2) + pow((uy - ty), 2))/2;}
float theta(int tx, int ty, int ux, int uy){return atan((uy - ty)/(ux - tx));}


int main(int argc, char **argv){

    // terminal input
    int opt;
    char* filename;
    int alpha = 0, percentage = 0, betha = 0, levelOneThreads = 0, levelTwoThreads = 0;


    // get options
    while ((opt = getopt(argc, argv, "f:a:r:b:u:d:")) != -1)
    {
        switch(opt){
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
            case '?':
                std::cerr << "Unknown option: " << char(optopt) << std::endl;
                return 1;
            default: 
                std::cerr << "Use: " << argv[0] << " -f <filename> -a <minimun alpha> -r <relative vote> -b <number of bethas> -u <number of threads 1> -d <number of threads 2>" << std::endl;
        }

        printf("%d", alpha);
    }

    // image variables
    fitsfile *fptr;
    int status;
    long fpixel = 1, naxis = 2, nelements, exposure;
    long naxes[2];

    status = 0;
    fits_open_file(&fptr, filename, READONLY, &status);
    fits_get_img_size(fptr, 2, naxes, &status);
    double *image = (double *) malloc(naxes[0]*naxes[1]*sizeof(double));

    fits_read_img(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1], NULL, image, NULL, &status);
    




    return 0;
}