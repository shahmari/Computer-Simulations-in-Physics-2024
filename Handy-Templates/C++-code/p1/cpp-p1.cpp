#include <iostream>
#include <math.h>
#include <string>
#include "./gnuplot-iostream.h"

using namespace std;


int main(){

    // We initilize arrays for our variables
    double x_array[100];
    double y_array[100];
    double z_array[100][100];

    // Using gnuplot for plotting
    Gnuplot gp1;
    
    // First we plot y = Sin(2x)
    // Calculating the function for each point
    for (int i=0; i<100; i++){
        x_array[i] = double(i) / 20.0;
        y_array[i] = sin(2.0*x_array[i]);

    }

    // Ploting and exportnig
    gp1 << "set xrange[0:5]\n set yrange[-1.5:1.5] \n";
    gp1 << "set size square \n";
    gp1 << "set title " << "\"y = Sin(2x) \"" <<"\n";
    gp1 << "set output " << "\"sin-plot.png\"" << "\n"; // Comment this line and the next for plot only mode
    gp1 << "set terminal pngcairo size 8000, 8000 fontscale 10 linewidth 10 pointscale 1\n";
    gp1 << "plot '-' with lines title \"\" \n";

    for (int i=0; i<100; i++){
        gp1 << to_string(x_array[i]) << " " << to_string(y_array[i]) << "\n";
    }

    // gp << "clear";

    // Now we calculate and plot z = Tanh(x*y)
    for (int i=0; i<100; i++){
        x_array[i] = double(i) / 20.0 - 2.5;
        for (int j=0; j<100; j++){
                y_array[j] = double(j) / 20.0 - 2.5;
                z_array[i][j] = tanh(x_array[i]*y_array[j]);
        }

    }



    // Ploting and exportnig
    Gnuplot gp2;
    gp2 << "set xrange[-2.5:2.5]\n set yrange[-2.5:2.5] \n";
    gp2 << "set size square \n";
    gp2 << "set title " << "\"z = Tanh(x*y) \"" <<"\n";
    gp2 << "set output " << "\"surface-plot.png\"" << "\n"; // Comment this line and the next for plot only mode
    gp2 << "set terminal pngcairo size 8000, 8000 fontscale 10 linewidth 10 pointscale 1\n";
    gp2 << "set hidden3d \n ";
    gp2 << "splot '-' with lines lc rgb \'#b90046\' lw 1.5 title \"\"\n";

    for (int i=0; i<100; i++){
        for (int j=0; j<100; j++){
            gp2 << to_string(x_array[i]) << " " << to_string(y_array[j]) << " " << to_string(z_array[i][j]) << "\n";
        }
        gp2 << "\n";
    }

    

    return 0;
}