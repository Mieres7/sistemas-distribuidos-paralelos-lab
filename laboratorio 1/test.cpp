#include <iostream>
#include <fitsio.h>

void check_status(int status) {
    if (status) {
        fits_report_error(stderr, status); // Informar del error de CFITSIO
        exit(status); // Salir si hay un error
    }
}

int main() {
    fitsfile *fptr;  // Puntero al archivo FITS
    int status = 0;  // Inicializar el estado como 0 (sin errores)
    long naxes[2] = {300, 300};  // Dimensiones de la imagen (300x300 píxeles)
    int bitpix = SHORT_IMG;  // Definir el tipo de datos de la imagen (16 bits)

    // Crear un nuevo archivo FITS
    fits_create_file(&fptr, "!testfile.fits", &status);
    check_status(status);

    // Crear la estructura básica de la imagen (una sola extensión)
    fits_create_img(fptr, bitpix, 2, naxes, &status);
    check_status(status);

    // Cerrar el archivo FITS
    fits_close_file(fptr, &status);
    check_status(status);

    std::cout << "Archivo FITS 'testfile.fits' creado con éxito." << std::endl;

    return 0;
}
