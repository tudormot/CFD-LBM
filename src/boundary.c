#include "boundary.h"
#include "computeCellValues.h"
#include "helper.h"
#include "LBDefinitions.h"

int is_valid(int x, int y, int z, dimensions dim){
    return (x >= 0) && (x <= dim.xlen + 1) && (y >= 0) && (y <= dim.ylen + 1) && (z >= 0) && (z <= dim.zlen + 1);
}

void treatBoundary(double *collideField_f, double *collideField_g, unsigned int* flagField, const double * const wallVelocity, dimensions dim){
    
    int x, y, z, xinv, yinv, zinv, idx, idxinv, f;
    int Q = 19;
    int xl = dim.xlen + 2;
    int xlyl = (xl)*(dim.ylen+2);
    double density, inner, weight_sum, T_d, T_local, Gc;
    
    // Velocity Boundary Conditions Looping in the faces
    // z-faces
    for(z = 0; z <= dim.zlen + 1; z+= (dim.zlen + 1) ){
        for(y = 0; y <= dim.ylen + 1; y++) {
            for(x = 0; x <= dim.xlen + 1; x++) {
                
                idx = z*xlyl + y*(xl) + x;
                
                f = flagField[idx];
                
                for(int i = 0; i < Q; i++) {
                    //Define inverse block
                    xinv = x + LATTICEVELOCITIES[i][0];
                    yinv = y + LATTICEVELOCITIES[i][1];
                    zinv = z + LATTICEVELOCITIES[i][2];
                    idxinv = zinv*xlyl + yinv*(xl) + xinv;
                    
                    if(is_valid(xinv, yinv, zinv, dim) && is_fluid(flagField[idxinv])){ //If the inverse cell is in the boundary or the domain
                        if (is_noslip(f)){
                            collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i];
                        }
                        else if (is_inflow(f)){
                            computeDensity(&collideField_f[Q*idxinv], &density); //Density of the inverse cell
                            inner = LATTICEVELOCITIES[i][0]*wallVelocity[0] + LATTICEVELOCITIES[i][1]*wallVelocity[1] 
                            + LATTICEVELOCITIES[i][2]*wallVelocity[2]; // dot_product(ci, Uwall)
                            
                            collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i] + 2.0*LATTICEWEIGHTS[i]*density*inner/C_S/C_S;
                            
                        }
                        else if (is_freeslip(f)){
                        collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i];
                        }
                    }
                }
                
                
                
            }
            
        }
    }
    
    // y-faces
    for(z = 0; z <= dim.zlen + 1; z++ ){
        for(y = 0; y <= dim.ylen + 1; y+= (dim.ylen + 1)) {
            for(x = 0; x <= dim.xlen + 1; x++) {
                
                idx = z*xlyl + y*(xl) + x;
                
                f = flagField[idx];
                
                for(int i = 0; i < Q; i++) {
                    //Define inverse block
                    xinv = x + LATTICEVELOCITIES[i][0];
                    yinv = y + LATTICEVELOCITIES[i][1];
                    zinv = z + LATTICEVELOCITIES[i][2];
                    idxinv = zinv*xlyl + yinv*(xl) + xinv;
                    
                    if(is_valid(xinv, yinv, zinv, dim) && is_fluid(flagField[idxinv])){ //If the inverse cell is in the boundary or the domain
                        if (is_noslip(f)){
                            collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i];
                        }
                        else if (is_inflow(f)){
                            computeDensity(&collideField_f[Q*idxinv], &density); //Density of the inverse cell
                            inner = LATTICEVELOCITIES[i][0]*wallVelocity[0] + LATTICEVELOCITIES[i][1]*wallVelocity[1] 
                            + LATTICEVELOCITIES[i][2]*wallVelocity[2]; // dot_product(ci, Uwall)
                            
                            collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i] + 2.0*LATTICEWEIGHTS[i]*density*inner/C_S/C_S;
                            
                        }
                        else if (is_freeslip(f)){
                        collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i];
                        }
                    }
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
                
                for(int i = 0; i < Q; i++) {
                    //Define inverse block
                    xinv = x + LATTICEVELOCITIES[i][0];
                    yinv = y + LATTICEVELOCITIES[i][1];
                    zinv = z + LATTICEVELOCITIES[i][2];
                    idxinv = zinv*xlyl + yinv*(xl) + xinv;
                    
                    if(is_valid(xinv, yinv, zinv, dim) && is_fluid(flagField[idxinv])){ //If the inverse cell is in the boundary or the domain
                        if (is_noslip(f)){
                            collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i];
                        }
                        else if (is_inflow(f)){
                            computeDensity(&collideField_f[Q*idxinv], &density); //Density of the inverse cell
                            inner = LATTICEVELOCITIES[i][0]*wallVelocity[0] + LATTICEVELOCITIES[i][1]*wallVelocity[1] 
                            + LATTICEVELOCITIES[i][2]*wallVelocity[2]; // dot_product(ci, Uwall)
                            
                            collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i] + 2.0*LATTICEWEIGHTS[i]*density*inner/C_S/C_S;
                            
                        }
                        else if (is_freeslip(f)){
                        collideField_f[Q*idx+i] = collideField_f[Q*idxinv+18-i];
                        }
                    }
                }
                
                
                
            }
            
        }
    }
    
    
    // Temperature Boundary Conditions
    for(z = 0; z <= dim.zlen + 1; z++ ){
        for(y = 0; y <= dim.ylen + 1; y++) {
            for(x = 0; x <= dim.xlen + 1; x++) {
                idx = z*xlyl + y*(xl) + x;
                
                f = flagField[idx];
                if (is_dirichl(f)){
                    get_sum_of_weights(x, y, z, dim, flagField, &weight_sum);               //Denominator of Gc expression
                    computeDensity(&collideField_f[Q*idx], &density);       //Density of current cell
                    get_Boundary_Temperature(x, y, z, &T_d);                //Dirichlet Temperature of current cell
                    computeTemperature(&collideField_g[Q*idx], &density, &T_local);   //Actual Temperature of current cell
                    
                    Gc = density*(T_d-T_local)/weight_sum;
                    
                    for(int i = 0; i < Q; i++) {
                        //Define inverse block
                        xinv = x + LATTICEVELOCITIES[i][0];
                        yinv = y + LATTICEVELOCITIES[i][1];
                        zinv = z + LATTICEVELOCITIES[i][2];
                        
                        idxinv = zinv*xlyl + yinv*(xl) + xinv;
                        if (is_fluid(flagField[idxinv])){          //If the inverse cell is not fluid
                            collideField_g[i] = LATTICEWEIGHTS[i]*Gc;
                        }
                    }
                }
                    
            }
        }
    }
    
}

void get_sum_of_weights(int x, int y, int z, dimensions dim, unsigned int* flagField, double* weight_sum){
    int xl = dim.xlen + 2;
    int xlyl = (xl)*(dim.ylen+2);
    int Q = 19;
    (*weight_sum) = 0.0;
    int xinv, yinv, zinv, idxinv;
    for(int i = 0; i < Q; i++) {
        //Define inverse block
        xinv = x + LATTICEVELOCITIES[i][0];
        yinv = y + LATTICEVELOCITIES[i][1];
        zinv = z + LATTICEVELOCITIES[i][2];
        
        idxinv = zinv*xlyl + yinv*(xl) + xinv;
        if (is_fluid(flagField[idxinv])){           //If the inverse cell is not fluid
            (*weight_sum) += LATTICEWEIGHTS[i];
        }
        
    }
}

void get_Boundary_Temperature(int x, int y, int z, double* T_d){
    
    if (x==1) (*T_d) = 0.0;
    else (*T_d) = 0.0;
    
}