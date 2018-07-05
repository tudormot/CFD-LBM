#include "boundary.h"
#include "computeCellValues.h"
#include "helper.h"

int is_valid(int x, int y, int z, dimensions dim){
    return (x >= 0) && (x <= dim.xlen + 1) && (y >= 0) && (y <= dim.ylen + 1) && (z >= 0) && (z <= dim.zlen + 1);
}

void treatBoundary(double *collideField_f,double *collideField_g, unsigned int* flagField, const double * const wallVelocity, dimensions dim){
    
    int x, y, z, xinv, yinv, zinv, idx, idxinv, f;
    int Q = 19;
    int xl = dim.xlen + 2;
    int xlyl = (xl)*(dim.ylen+2);
    double density, inner;
    
    // z-faces
    for(z = 0; z <= dim.zlen + 1; z+= (dim.zlen + 1) ){
        for(y = 0; y <= dim.ylen + 1; y++) {
            for(x = 0; x <= dim.xlen + 1; x++) {
                
                idx = z*xlyl + y*(xl) + x;
                
                f = flagField[idx];
                switch (f){
                    
                    case 1: // No Slip
                    for(int i = 0; i < Q; i++) {
                        //Define inverse block
                        xinv = x + LATTICEVELOCITIES[i][0];
                        yinv = y + LATTICEVELOCITIES[i][1];
                        zinv = z + LATTICEVELOCITIES[i][2];
                        
                        if(is_valid(xinv, yinv, zinv, dim)){ //If the inverse cell is in the boundary or the domain
                            idxinv = zinv*xlyl + yinv*(xl) + xinv;
                            if (is_fluid(flagField[idxinv])){          //If the inverse cell is fluid
                                collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i];
                            }
                        }
                    }
                    break;
                    
                    case 2: // Moving Wall
                    for(int i = 0; i < Q; i++) {
                        //Define inverse block
                        xinv = x + LATTICEVELOCITIES[i][0];
                        yinv = y + LATTICEVELOCITIES[i][1];
                        zinv = z + LATTICEVELOCITIES[i][2];
                        
                        if(is_valid(xinv, yinv, zinv, dim)){ //If the inverse cell is in the boundary or the domain
                            idxinv = zinv*xlyl + yinv*(xl) + xinv;
                            if (flagField[idxinv]==0){          //If the inverse cell is fluid
                                
                                computeDensity(&collideField_f[Q*idxinv], &density); //Density of the inverse cell
                                inner = LATTICEVELOCITIES[i][0]*wallVelocity[0] + LATTICEVELOCITIES[i][1]*wallVelocity[1] 
                                        + LATTICEVELOCITIES[i][2]*wallVelocity[2]; // dot_product(ci, Uwall)
                                
                                collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i] + 2.0*LATTICEWEIGHTS[i]*density*inner/C_S/C_S;
                            }
                        }
                    }
                    break;
                    
                }
                
            }
        }
    }
    
    // y-faces
    for(z = 0; z <= dim.zlen + 1; z++ ){
        for(y = 0; y <= dim.ylen + 1; y+= (dim.ylen+ 1)) {
            for(x = 0; x <= dim.xlen + 1; x++) {
                
                idx = z*xlyl + y*(xl) + x;
                
                f = flagField[idx];
                switch (f){
                    
                    case 1: // No Slip
                    for(int i = 0; i < Q; i++) {
                        //Define inverse block
                        xinv = x + LATTICEVELOCITIES[i][0];
                        yinv = y + LATTICEVELOCITIES[i][1];
                        zinv = z + LATTICEVELOCITIES[i][2];
                        
                        if(is_valid(xinv, yinv, zinv, dim)){ //If the inverse cell is in the boundary or the domain
                            idxinv = zinv*xlyl + yinv*(xl) + xinv;
                            if (is_fluid(flagField[idxinv])){          //If the inverse cell is fluid
                                collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i];
                            }
                        }
                    }
                    break;
                    
                    case 2: // Moving Wall
                    for(int i = 0; i < Q; i++) {
                        //Define inverse block
                        xinv = x + LATTICEVELOCITIES[i][0];
                        yinv = y + LATTICEVELOCITIES[i][1];
                        zinv = z + LATTICEVELOCITIES[i][2];
                        
                        if(is_valid(xinv, yinv, zinv, dim)){ //If the inverse cell is in the boundary or the domain
                            idxinv = zinv*xlyl + yinv*(xl) + xinv;
                            if (is_fluid(flagField[idxinv])){          //If the inverse cell is fluid
                                
                                computeDensity(&collideField_f[Q*idxinv], &density); //Density of the inverse cell
                                inner = LATTICEVELOCITIES[i][0]*wallVelocity[0] + LATTICEVELOCITIES[i][1]*wallVelocity[1] 
                                        + LATTICEVELOCITIES[i][2]*wallVelocity[2]; // dot_product(ci, Uwall)
                                
                                collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i] + 2.0*LATTICEWEIGHTS[i]*density*inner/C_S/C_S;
                            }
                        }
                    }
                    break;
                    
                }
                
            }
        }
    }
    
    // x-faces
    for(z = 0; z <= dim.zlen + 1; z++ ){
        for(y = 0; y <= dim.ylen + 1; y++) {
            for(x = 0; x <= dim.xlen + 1; x+= (dim.xlen + 1)) {
                
                idx = z*xlyl + y*(xl) + x;
                
                f = flagField[idx];
                switch (f){
                    
                    case 1: // No Slip
                    for(int i = 0; i < Q; i++) {
                        //Define inverse block
                        xinv = x + LATTICEVELOCITIES[i][0];
                        yinv = y + LATTICEVELOCITIES[i][1];
                        zinv = z + LATTICEVELOCITIES[i][2];
                        
                        if(is_valid(xinv, yinv, zinv, dim)){ //If the inverse cell is in the boundary or the domain
                            idxinv = zinv*xlyl + yinv*(xl) + xinv;
                            if (is_fluid(flagField[idxinv])){          //If the inverse cell is fluid
                                collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i];
                            }
                        }
                    }
                    break;
                    
                    case 2: // Moving Wall
                    for(int i = 0; i < Q; i++) {
                        //Define inverse block
                        xinv = x + LATTICEVELOCITIES[i][0];
                        yinv = y + LATTICEVELOCITIES[i][1];
                        zinv = z + LATTICEVELOCITIES[i][2];
                        
                        if(is_valid(xinv, yinv, zinv, dim)){ //If the inverse cell is in the boundary or the domain
                            idxinv = zinv*xlyl + yinv*(xl) + xinv;
                            if (is_fluid(flagField[idxinv])){          //If the inverse cell is fluid
                                
                                computeDensity(&collideField_f[Q*idxinv], &density); //Density of the inverse cell
                                inner = LATTICEVELOCITIES[i][0]*wallVelocity[0] + LATTICEVELOCITIES[i][1]*wallVelocity[1] 
                                        + LATTICEVELOCITIES[i][2]*wallVelocity[2]; // dot_product(ci, Uwall)
                                
                                collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i] + 2.0*LATTICEWEIGHTS[i]*density*inner/C_S/C_S;
                            }
                        }
                    }
                    break;
                    
                }
                
            }
        }
    }
    
}
