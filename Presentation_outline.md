PRESENTATION PLAN

1. Introduction - What are we going to talk about?
	- We did LBM driven cavity and heat transport in 3D
	- Aim of the presentaion is to compare NSE with LBM from different perspectives
    - Show application from industry (nice pictures Exa etc.,
2. Main Part
	2.1 Introduction to LBM - theory (including stability etc. ) and comparison to NSE in summary. (Only theory here).
		2.1.1 Boltzmann equation
        2.1.3 Relationship to macroscopic parameters of interest
        2.1.2 LatticeBoltzmann method(after the discretization)
        2.3.x Linearization of colision operator (BGK)
        2.3.x Theory related to heat:
            - very similar to the velocity calculations                
            - energy distribution
            - convection required gravity: we had to use the Boussinesq approx, which gets added to the velocity field
            - bla bla whatever we find interesting in the journal when rereading    
            - as for the velocity case, we have 2 possible boundary types: Dorichlet and Neumann, but we will get into detail into the implementation part of the presentation
 
    2.2 Implementation - #show main in pseudo code, talk about functions like doColiision, treatBoundary etc., mention flags (binary stuff).
        -general information on the memory structures used: with emphasis on streaming field and collision field
        -splitting of the lattice boltzmann equations with bgk collision into 2 steps for ease and performance of computation
        -boundary tratement:
            -velocity: Dirichlet, Neumann
            -temperature:   Dirichlet, Neumann
        -show pseudocode of main for velocity without temperature: make the point that there is no solving step, everything is explicit( we'l come back to that in the results phase)
        -show pesudocode of main with both velocity and temperature: explain what was added
    2.3 Results
	    2.2.1 Results for Driven Cavity and comparison
		    - methodology (1-slice, free slip)
            - show pictures of both LBM and NS
            - explain why the pictures are different: we know that they are pictures, they have to do with the very different parameter settings in LBM and NS, and this is actually one of the disadvantages of LBM, the fact that the parameters are hard to be chosen..
            - show runtime difference of LBM and NS
            - memory consumption when running the two cases
	    2.2..2 Results for Heat Transport and comparison
		    - do the same thing as above
    2.4 Comparison more detailed:
        2.1.4 where is NSE typicaly used vs LBM - 
        -difference in theory: LBM has roots in statistics whereas NS is different bla bla..
        - Advantages and disadvantages of LBM:
            -LBM does not contain a solving step. Hence , it is a purely explicit method, stabiity is controled via mesh density
            -LBM much easier/cleaner to implement (once you understand the theory): here show 2 pictures together (of pseudocode of NS and pseud. of LBM) 
            -LBM has this problem of parameters not having much physical significance
            -LBM is a very local method, hence very easy to implement
            -NS has more adaptability: it is trivial to refine the mesh only locally , withou much changes to the algorithm. At leas in our simple LBM implementation, this would not be possible, we would have to refine the whole grid
            -LBM works well for weakly compressible flows. Not usable for incompressible flows, as everything relies of density distribution. Does not work well for lows with high Mach number. 
            -Does not function well for turbulence (probably becuase of lack of implicit step (doublecheck :) )) 
            
	

    
3. Summary - most important things that we want people to remember, not too long. E.g. 
    todo: LATER


	

	

