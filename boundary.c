#include "boundary.h"


void setZone2DBoundaryFlags(struct gridZone2D zone)
{
    /* Zone is at the left end of a tile */
    if (zone.iTile == 0)
    {
        /* Is the zone at the left physical boundary? */
        if (zone.i == 0)
        {   
            /* If yes, then apply the physical boundary flag */

             zone.leftEdge = PHYSICAL_BOUNDARY_LEFT_EDGE; 
        }
        else
        {
            /* Zone is at the tile/MPI boundary but not at the physical boundary */

            zone.leftEdge = GET_DATA_FROM_GLOBAL_ARRAY;
        }
    }

    /* Zone is at the right end of a tile */
    if (zone.iTile == TILE_SIZE_X1-1)
    {
        /* Is the zone at the right physical boundary? */
        if (zone.i == N1-1)
        {
            /* If yes, then apply the physical boundary condition flag */

            zone.rightEdge = PHYSICAL_BOUNDARY_RIGHT_EDGE;
        }
        else
        {   
            /* Zone is at the tile/MPI boundary but not at the physical boundary */

            zone.rightEdge = GET_DATA_FROM_GLOBAL_ARRAY;
        }

    }

    /* Zone is at the top end of a tile */
    if (zone.jTile == 0)
    {
        /* Is the zone at the top physical boundary? */
        if (zone.j == 0)
        {
            /* If yes, the apply the physical boundary flag */

            zone.topEdge = PHYSICAL_BOUNDARY_TOP_EDGE;
        }
        else
        {
            /* Zone is at the tile/MPI boundary but not at the physical boundary */

            zone.topEdge = GET_DATA_FROM_GLOBAL_ARRAY;
        }
    }

    /* Zone is at the bottom end of a tile */
    if (zone.jTile == TILE_SIZE_X2-1)
    {
        /* Is the zone at the bottom physical boundary? */
        if (zone.j == N2-1)
        {
            /* If yes, then apply the physical boundary flag */

            zone.bottomEdge = PHYSICAL_BOUNDARY_BOTTOM_EDGE;
        }
        else
        {
            /* Zone is at the tile/MPI boundary but not at the physical boundary */

            zone.bottomEdge = GET_DATA_FROM_GLOBAL_ARRAY;
        }
    }

}


void applyTileBoundaryConditions(const struct gridZone2D zone,
                                 const REAL restrict *globalPrimArray,
                                 REAL restrict *tile)
{

    /* LEFT EDGE */
    if (zone.leftEdge == GET_DATA_FROM_GLOBAL_ARRAY)
    {
        for (int var=0; var<DOF: var++)
        {
            for (int iGhost=-NG; iGhost<0; iGhost++)
            {
                /* Get data from globalPrimArray */
                tile[INDEX_TILE(iGhost, zone.jTile, var)] = 
                globalPrimArray[INDEX_GLOBAL(zone.i+iGhost, zone.j, var)];
            }
        }
    }
    else if (zone.leftEdge == OUTFLOW)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int iGhost=-NG; iGhost<0; iGhost++)
            {
                /* Copy the value of the variables at the left edge of the tile.
                 * Note that the left edge of the tile when zone.leftEdge is not
                 * TILE_BOUNDARY coincides with the left boundary.*/
                tile[INDEX_TILE(iGhost, zone.jTile, var)] =
                tile[INDEX_TILE(0, zone.jTile, var)];
            }
        }
    }
    else if (zone.leftEdge == MIRROR)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int iGhost=-NG; iGhost<0; iGhost++)
            {
                /* Mirror the left edge of the tile */
                tile[INDEX_TILE(iGhost, zone.jTile, var)] = 
                tile[INDEX_TILE(-iGhost-1, zone.jTile, var)];
            }
        }
    }
    else if (zone.leftEdge == PERIODIC_SINGLE_MPI_NODE)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int iGhost=-NG; iGhost<0; iGhost++)
            {
                /* Copy data from the right edge of the global domain */
                tile[INDEX_TILE(iGhost, zone.jTile, var)] = 
                globalPrimArray[INDEX_GLOBAL(N1+iGhost, zone.j, var)];
            }
        }
    }
    /* END OF LEFT EDGE */

    /* RIGHT EDGE */
    if (zone.rightEdge == GET_DATA_FROM_GLOBAL_ARRAY)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int iGhost=0; iGhost<NG; iGhost++)
            {
                tile[INDEX_TILE(TILE_SIZE_X1+iGhost, zone.jTile, var)] = 
                globalPrimArray[INDEX_GLOBAL(zone.i+1+iGhost, zone.j, var)];
            }
        }
    }
    else if (zone.rightEdge == OUTFLOW)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int iGhost=0; iGhost<NG; iGhost++)
            {
                tile[INDEX_TILE(TILE_SIZE_X1+iGhost, zone.jTile, var)] =
                tile[INDEX_TILE(TILE_SIZE_X1-1, zone.jTile, var)];
            }
        }
    }
    else if (zone.rightEdge == MIRROR)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int iGhost=0; iGhost<NG; iGhost++)
            {
                tile[INDEX_TILE(TILE_SIZE_X1+iGhost, zone.jTile, var)] =
                tile[INDEX_TILE(TILE_SIZE_X1-1-iGhost, zone.jTile, var)];
            }
        }

    }
    else if (zone.rightEdge == PERIODIC_SINGLE_MPI_NODE)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int iGhost=0; iGhost<NG; iGhost++)
            {
                tile[INDEX_TILE(TILE_SIZE_X1+iGhost, zone.jTile, var)] = 
                globalPrimArray[INDEX_GLOBAL(iGhost, zone.j, var)];
            }
        }
    }
    /* END OF RIGHT EDGE */

    /* TOP EDGE */
    if (zone.topEdge == GET_DATA_FROM_GLOBAL_ARRAY)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int jGhost=-NG; jGhost<0; jGhost++)
            {
                tile[INDEX_TILE(zone.iTile, jGhost, var)] = 
                globalPrimArray[INDEX_GLOBAL(zone.i, zone.j+jGhost, var)];
            }
        }
    }
    else if (zone.topEdge == OUTFLOW)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int jGhost=-NG; jGhost<0; jGhost++)
            {
                tile[INDEX_TILE(zone.iTile, jGhost, var)] =
                tile[INDEX_TILE(zone.iTile, 0, var)];
            }
        }
    }
    else if (zone.topEdge == MIRROR)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int jGhost=-NG; jGhost<0; jGhost++)
            {
                tile[INDEX_TILE(zone.iTile, jGhost, var)] = 
                tile[INDEX_TILE(zone.iTile, -jGhost-1, var)];
            }
        } 
    }
    else if (zone.topEdge == PERIODIC_SINGLE_MPI_NODE)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int jGhost=-NG; jGhost<0; jGhost++)
            {
                tile[INDEX_TILE(zone.iTile, jGhost, var)] = 
                globalPrimArray[INDEX_GLOBAL(zone.i, N2+jGhost, var)];
            }
        }
    }
    /* END OF TOP EDGE */

    /* BOTTOM EDGE */
    if (zone.bottomEdge == GET_DATA_FROM_GLOBAL_ARRAY)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int jGhost=0; jGhost<NG; jGhost++)
            {
                tile[INDEX_TILE(zone.iTile, TILE_SIZE_X2+jGhost, var)] = 
                globalPrimArray[INDEX_GLOBAL(zone.i, zone.j+1+jGhost, var)];
            }
        }
    }
    else if (zone.bottomEdge == OUTFLOW)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int jGhost=0; jGhost<NG; jGhost++)
            {
                tile[INDEX_TILE(zone.iTile, TILE_SIZE_X2+jGhost, var)] = 
                tile[INDEX_TILE(zone.iTile, TILE_SIZE_X2-1, var)];
            }
        }
    }
    else if (zone.bottomEdge == MIRROR)
    {
        for (int var=0; var<DOF; var++)
        {
            for (int jGhost=0; jGhost<NG; jGhost++)
            {
                tile[INDEX_TILE(zone.iTile, TILE_SIZE_X2+jGhost, var)] = 
                tile[INDEX_TILE(zone.iTile, TILE_SIZE_X2-1-jGhost, var)];
            }
        }
    }
    else if (zone.bottomEdge == PERIODIC_SINGLE_MPI_NODE)
    {   
        for (int var=0; var<DOF; var++)
        {
            for (int jGhost=0; jGhost<NG; jGhost++)
            {
                tile[INDEX_TILE(zone.iTile, TILE_SIZE_X2+jGhost, var)] = 
                globalPrimArray[INDEX_GLOBAL(zone.i, jGhost, var)];
            }
        }
    }
    /* END OF BOTTOM EDGE */

}

