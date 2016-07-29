#ifndef GRIM_RECONSTRUCT_H_
#define GRIM_RECONSTRUCT_H_

#include "../grid/grid.hpp"

/* Reconstruction routines */
namespace reconstruction
{
  array minmod(array &x, array &y, array &z,
               int &numReads, int &numWrites
              );
  
  array slopeMM(const int dir,const double dX,
		            const array& in,
                int &numReads,
                int &numWrites
               );

  array slopeWENO5(const int dir,const double dX,
		               const array& in,
                   int &numReads,
                   int &numWrites
                  );

  array slopePPM(const int dir,const double dX,
		               const array& in,
                   int &numReads,
                   int &numWrites
                  );

  array slope(const int dir,const double dX,
		          const array& in,
              int &numReads,
              int &numWrites
             );

  void reconstructMM(const grid &prim,
                     const int dir,
	              		 grid &primLeft,
                     grid &primRight,
                     int &numReads,
                     int &numWrites
                    );

  void reconstructWENO5(const grid &prim,
                  			const int dir,
                  			grid &primLeft,
                  			grid &primRight,
                        int &numReads,
                        int &numWrites
                  		 );

  void reconstructPPM(const grid &prim,
                  			const int dir,
                  			grid &primLeft,
                  			grid &primRight,
                        int &numReads,
                        int &numWrites
                  		 );

  void reconstruct(const grid &prim,
                   const int dir,
                   grid &primLeft,
                   grid &primRight,
                   int &numReads,
                   int &numWrites
                  );
}

#endif /* GRIM_RECONSTRUCT_H_ */
