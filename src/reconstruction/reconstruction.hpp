#ifndef GRIM_RECONSTRUCT_H_
#define GRIM_RECONSTRUCT_H_

#include "../grid/grid.hpp"

/* Reconstruction routines */
namespace reconstruction
{
  array minmod(array &x, array &y, array &z);
  
  array slopeMM(const int dir,const double dX,
		            const array& in
               );

  void reconstructMM(const grid &prim,
                     const int dir,
	              		 grid &primLeft,
                     grid &primRight
                    );

  void reconstructWENO5(const grid &prim,
                  			const int dir,
                  			grid &primLeft,
                  			grid &primRight
                  		 );

  void reconstruct(const grid &prim,
                   const int dir,
                   grid &primLeft,
                   grid &primRight
                  );
}

#endif /* GRIM_RECONSTRUCT_H_ */
