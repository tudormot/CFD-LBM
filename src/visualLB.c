#include "visualLB.h"
#include "helper.h"
#include "LBDefinitions.h"
#include <math.h>
#include "computeCellValues.h"

void writeVtkOutput(const double * const collideField, const unsigned int * const flagField,
                    const char* filename, int t, int xlength, double* vel){
    
    int x, y, z; //Indices for x, y, z and velocity respectively
   // int i; // commented, because unused
   // int Q = 19;
    //double density;
   // double velocity[3];
  //  const double* currentCell;
    int xl2 = (xlength+2)*(xlength+2);
    int idx;

    // Open File
    char szFileName[256];
    FILE *fp=NULL;
    sprintf( szFileName, "Output/Out_%i.vtk", t);
    fp = fopen( szFileName, "w");
    if( fp == NULL ){
        char szBuff[256];
        sprintf( szBuff, "Failed to open %s", szFileName);
        ERROR( szBuff );
        return;
    }
    
    // Write Header and Coordinates
    write_vtkHeader(fp, xlength);
    write_vtkPointCoordinates(fp, xlength);
    
    // Write Velocities
    fprintf(fp, "\nPOINT_DATA %i \n", xlength*xlength*xlength);
    fprintf(fp,"\n");
    fprintf(fp, "VECTORS velocity float\n");
    
    for(z = 1; z <= xlength; z++){
        for(y = 1; y <= xlength; y++) {
            for(x = 1; x <= xlength; x++) {
                idx = (z*xl2 + y*(xlength+2) + x);

                    fprintf(fp, "%f %f %f\n", vel[3*idx+0], vel[3*idx+1], vel[3*idx+2] );

            }
        }
    }
    
    // Close File
    if( fclose(fp) )
    {
        char szBuff[80];
        sprintf( szBuff, "Failed to close %s", szFileName );
        ERROR( szBuff );
    }
    
}
        
void write_vtkHeader( FILE *fp, int xlength) {
    if( fp == NULL )		       
    {
        char szBuff[80];
        sprintf( szBuff, "Null pointer in write_vtkHeader" );
        ERROR( szBuff );
        return;
    }
    
    fprintf(fp,"# vtk DataFile Version 2.0\n");
    fprintf(fp,"generated by CFD-lab course output\n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"\n");	
    fprintf(fp,"DATASET STRUCTURED_GRID\n");
    fprintf(fp,"DIMENSIONS  %i %i %i \n", xlength, xlength, xlength);
    fprintf(fp,"POINTS %i int\n", (xlength)*(xlength)*(xlength) );
    fprintf(fp,"\n");
}

void write_vtkPointCoordinates( FILE *fp, int xlength) {

    int x, y, z;
    for(z = 1; z <= xlength; z++){
        for(y = 1; y <= xlength; y++) {
            for(x = 1; x <= xlength; x++) {
                fprintf(fp, "%i %i %i\n", x, y, z );
            }
        }
    }
}
