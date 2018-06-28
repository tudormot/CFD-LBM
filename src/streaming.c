#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){

	//TODO: if code is too slow this function is a good candidate for changing, as it could be rewritten without the if..
	int x_dest, y_dest, z_dest; //coordinates of destination cell
	for(int z = 1; z <= xlength; z++ ){
		for(int y = 1; y <= xlength; y++) {
			for(int x = 1; x <= xlength ; x++) {
				for(int l = 0;l<NO_OF_LATTICE_DIMENSIONS;l++)
				{
					/*determine which is the destination cell based on current cell and lattice velocities*/
					x_dest = x + LATTICEVELOCITIES[l][0];
					y_dest = y + LATTICEVELOCITIES[l][1];
					z_dest = z + LATTICEVELOCITIES[l][2];

					/*check if destinaton is a valid destination first*/
					if( (x_dest >=0) && (y_dest >=0) && (z_dest >=0) && (x_dest <= xlength+1 ) && (y_dest <= xlength+1 ) && (z_dest <= xlength+1 ))
					{
						/*now perform the streaming step, from source to destination*/
						streamField[NO_OF_LATTICE_DIMENSIONS*(z_dest*xlength*xlength + y_dest*xlength + x_dest) + l]=\
								collideField[NO_OF_LATTICE_DIMENSIONS*(z*xlength*xlength + y*xlength + x) + l];
					}
				}

			}
		}
	}
}
