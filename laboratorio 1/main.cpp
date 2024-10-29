#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cmath> 
#include <omp.h>
#include <fitsio.h>
#include <vector>

using namespace std;



struct Ellipse {
    int ox, oy;
    double  alpha ,betha, tetha;
};

// Data functions
int getOx(int tx, int ux) { return (tx + ux) / 2; }
int getOy(int ty, int uy) { return (ty + uy) / 2; }
double getAlpha(int tx, int ty, int ux, int uy) { return (sqrt(pow((ux - tx), 2) + pow((uy - ty), 2)) / 2); }
double getTheta(int tx, int ty, int ux, int uy) { return atan2((uy - ty),(ux - tx)); }
double getDelta(int kx, int ky, int ox, int oy) { return sqrt(pow(ky-oy,2) + pow(kx-ox,2));}
double getGamma(double theta, int kx, int ky, int ox, int oy){ return sin(theta) * (ky - oy) + cos(theta) * (kx - ox); }
double getBetha(double alpha, double gamma, double delta) 
{ 
    return sqrt((pow(alpha, 2) * pow(delta,2) - (pow(alpha, 2) * pow(gamma, 2)))/ (pow(alpha, 2) - pow(gamma, 2))) ; 
}



int main(int argc, char **argv) {

    // terminal input
    int opt;
    char* filename = nullptr;
    int minAlpha = 0, numberOfBethas = 0, levelOneThreads = 0, levelTwoThreads = 0;
    float relativeVote = 0;

    // get options
    while ((opt = getopt(argc, argv, "i:a:r:b:u:d:")) != -1) {
        switch (opt) {
            case 'i':
                filename = optarg;
                break;
            case 'a':
                minAlpha = std::atoi(optarg);
                break;
            case 'r':
                relativeVote = std::atof(optarg);
                break;
            case 'b':
                numberOfBethas = std::atoi(optarg);
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
                std::cerr << "Use: " << argv[0] << " -i <filename> -a <minimun alpha> -r <relative vote> -b <number of bethas> -u <number of threads 1> -d <number of threads 2>" << std::endl;
        }
    }

    if (filename == nullptr) {
        cerr << "Error: Debe proporcionar un nombre de archivo con -i." << endl;
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
    // Delta Beta for discretizacion
    double deltaBetha = static_cast<double>(naxes[0]) / (2 * numberOfBethas);
    

    // First parallel level - create edge list 
    omp_set_nested(1);
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
        #pragma omp parallel num_threads(levelTwoThreads) shared(pixels, ellipses, naxes, deltaBetha, minAlpha, relativeVote)
        {
            #pragma omp for 
            for(pair<int,int> pixel_t: pixels ){

                // #pragma omp for
                for(pair<int, int> pixel_u : pixels){

                    // array for betha votes
                    vector<int> votes = vector<int>(numberOfBethas, 0);

                    int oX = getOx(pixel_t.first, pixel_u.first);
                    int oY = getOy(pixel_t.second, pixel_u.second);
                    double alpha = getAlpha(pixel_t.first, pixel_t.second, pixel_u.first, pixel_u.second);
                    double theta = getTheta(pixel_t.first, pixel_t.second, pixel_u.first, pixel_u.second);

                    if(alpha <= minAlpha){continue;}

                    for(pair<int, int> pixel_k : pixels){
                        
                        if(pixel_k == pixel_t || pixel_k == pixel_u){ continue; }

                        double delta = getDelta(pixel_k.first, pixel_k.second, oX, oY);
                        if(delta > alpha){ continue; }

                        double gamma = getGamma(theta, pixel_k.first, pixel_k.second, oX, oY);
                        double betha = getBetha(alpha, gamma, delta);
                        
                        if(isnan(betha) || betha <= 0 || betha > alpha || alpha > (max(naxes[0], naxes[1])) / 2){ continue; }

                        int bethaIndex = betha / deltaBetha;

                        #pragma omp atomic
                        votes[bethaIndex]++; 
                        
                    }
                    
                    #pragma omp critical                 
                    for(size_t i = 0; i < votes.size(); i++){
                        double b = i * deltaBetha;
                        int circumference = M_PI * (3*(alpha + b) - sqrt((3*alpha + b) * (alpha + 3*b)));
                        
                        if(votes[i] > relativeVote * circumference){
                            Ellipse ep = Ellipse{oX, oY, alpha, static_cast<double>(votes[i]),  theta*(180/M_PI)};
                            ellipses.insert(ellipses.end(), ep);
                            printf("%d %d %f %f %f\n", oX, oY, alpha, i*deltaBetha, theta*(180/M_PI));
                        }
                    }
                    
                }
            }
        }
    }
    
    
    printf("ellipses found: %li", ellipses.size());

    // Close fits file
    if (fits_close_file(fptr, &status)) {
        fits_report_error(stderr, status);
        free(image);
        return status;
    }



    free(image);

    return 0;
}
