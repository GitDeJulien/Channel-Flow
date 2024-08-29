#include <stdio.h>
#include <math.h>

#define pi acos(-1.0)

int main() {
    // Define the parameters
    char model_name[] = "WRLES_Retau395";
    float xlen = 2*pi;
    float ylen = pi;
    float zlen = 2.0;
    float rho = 1.0;
    float uinf = 1.0;
    float re = 13800.0;
    float ret = 388.97;
    float dt = 0.01;
    float g_limit0 = 3;
    float g_limit1 = 40;
    int zplus[] = {5, 20, 40, 60, 80, 98, 151, 199, 251, 302};
    char in_path[] = "/media/julien/E91C-A3F7/channel_wrles_retau395/";
    char out_figure_path[] = "output/figures/channel_wrles_retau395/";
    char rans_path[] = "input/chan_rans_mean.dat";
    char out_data_path[] = "output/datas/channel_wrles_retau395/";

    // Open the file for writing
    FILE *file = fopen("parameters.dat", "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    // Write the parameters to the file
    fprintf(file, "model_name = %s\n", model_name);
    fprintf(file, "xlen = %.4f\n", xlen);
    fprintf(file, "ylen = %.4f\n", ylen);
    fprintf(file, "zlen = %.2f\n", zlen);
    fprintf(file, "rho = %.2f\n", rho);
    fprintf(file, "uinf = %.2f\n", uinf);
    fprintf(file, "re = %.4f\n", re);
    fprintf(file, "ret = %.2f\n", ret);
    fprintf(file, "dt = %.4f\n", dt);
    fprintf(file, "g_limit0 = %.2f\n", g_limit0);
    fprintf(file, "g_limit1 = %.2f\n", g_limit1);
    
    // Write the zplus array
    fprintf(file, "zplus = [");
    for (int i = 0; i < sizeof(zplus)/sizeof(zplus[0]); i++) {
        if (i != 0) {
            fprintf(file, ", ");
        }
        fprintf(file, "%d", zplus[i]);
    }
    fprintf(file, "]\n");

    // Write the strings
    fprintf(file, "in_path = %s\n", in_path);
    fprintf(file, "out_figure_path = %s\n", out_figure_path);
    fprintf(file, "rans_path = %s\n", rans_path);
    fprintf(file, "out_data_path = %s\n", out_data_path);

    // Close the file
    fclose(file);

    printf("Parameters written to parameters.dat\n");
    return 0;
}

// gcc -o write_parameters_wrles_395 write_parameters_wrles_395.c