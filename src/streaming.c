#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,unsigned int *flagField,dimensions dim){

	int x_dest, y_dest, z_dest; //coordinates of destination cell
	
	int xl = dim.xlen + 2;
	int xlyl = xl*(dim.ylen+2);
	
	for(int z = 0; z <= dim.zlen+1; z++ ){
		for(int y = 0; y <= dim.ylen+1; y++) {
			for(int x = 0; x <= dim.xlen+1 ; x++) {
				for(int l = 0;l<NO_OF_LATTICE_DIMENSIONS;l++)
				{
					/*determine which is the destination cell based on current cell and lattice velocities*/
					x_dest = x + LATTICEVELOCITIES[l][0];
					y_dest = y + LATTICEVELOCITIES[l][1];
					z_dest = z + LATTICEVELOCITIES[l][2];

					/*check if destinaton is a valid destination first*/
					if( (x_dest >=1) && (y_dest >=1) && (z_dest >=1) && (x_dest <= dim.xlen ) && (y_dest <= dim.ylen ) && (z_dest <= dim.zlen ))
					{
						/*now perform the streaming step, from source to destination*/
						streamField[NO_OF_LATTICE_DIMENSIONS*(z_dest*xlyl + y_dest*xl + x_dest) + l]= collideField[NO_OF_LATTICE_DIMENSIONS*(z*xlyl + y*xl + x) + l];
					}
				}

			}
		}
	}
}
