#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,unsigned int *flagField,int xlength){

	//TODO: if code is too slow this function is a good candidate for changing, as it could be rewritten without the if..
	int x_dest, y_dest, z_dest; //coordinates of destination cell
	
	int xl = xlength + 2;
	int xl2 = xl*xl;
	
	for(int z = 0; z <= xlength+1; z++ ){
		for(int y = 0; y <= xlength+1; y++) {
			for(int x = 0; x <= xlength+1 ; x++) {
				for(int l = 0;l<NO_OF_LATTICE_DIMENSIONS;l++)
				{
					/*determine which is the destination cell based on current cell and lattice velocities*/
					x_dest = x + LATTICEVELOCITIES[l][0];
					y_dest = y + LATTICEVELOCITIES[l][1];
					z_dest = z + LATTICEVELOCITIES[l][2];

					/*check if destinaton is a valid destination first*/
					if( (x_dest >=1) && (y_dest >=1) && (z_dest >=1) && (x_dest <= xlength ) && (y_dest <= xlength ) && (z_dest <= xlength ))
					{
						/*now perform the streaming step, from source to destination*/
						streamField[NO_OF_LATTICE_DIMENSIONS*(z_dest*xl2 + y_dest*xl + x_dest) + l]= collideField[NO_OF_LATTICE_DIMENSIONS*(z*xl2 + y*xl + x) + l];
					}
				}

			}
		}
	}
}
