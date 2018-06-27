#include "boundary.h"

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
    
    int x, y, z, xinv, yinv, zinv, idx, f;
    int Q = 19;
    int xlength2 = (xlength+2)*(xlength+2);
    
    for(z = 0; z <= xlength + 1; z+= (xlength + 1) ){
        for(y = 0; y <= xlength + 1; y+= (xlength + 1)) {
            for(x = 0; x <= xlength + 1; x+= (xlength + 1)) {
                
                idx = z*xlength2 + y*(xlength+2) + x;
                
                f = flagField[idx];
                switch (f){
                    
                    case 1: // No Slip
                    
                    for(i = 0; i < Q; i++) {
                        //Define inverse block
                        xinv = x + LATTICEVELOCITIES[i][0];
                        yinv = y + LATTICEVELOCITIES[i][1];
                        zinv = z + LATTICEVELOCITIES[i][2];
                        
                        collideField[Q*idx+i] = 
                        break;
                        
                        case 2: // Moving Wall
                        break;
                        
                    }
                    
                }}
            }
        }
        
        
        
    }
    
