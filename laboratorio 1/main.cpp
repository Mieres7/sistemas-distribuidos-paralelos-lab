#include <iostream>
#include <unistd.h>
#include <cstdlib>

#include <fitsio.h>



int main(int argc, char **argv){

    // terminal input
    int opt;
    std::string filename;
    int alpha = 0, percentage = 0, betha = 0, threads = 0;

    // get options
    while ((opt = getopt(argc, argv, "f:a:r:b:T:")) != -1)
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
            case 'T':
                threads = std::atoi(optarg);
                break;
            case '?':
                std::cerr << "Unknown option: " << char(optopt) << std::endl;
                return 1;
            default: 
                std::cerr << "Use: " << argv[0] << " -f <filename> -a <alpha> -r <percetage apply to circumsphere> -b <number of bethas> -T <number of threads>" << std::endl;
        }

        printf("%d", alpha);


    }
    


    return 0;
}