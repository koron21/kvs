#include "BoundingBox.h"
#include <cmath>

namespace kvs
{
    BoundingBox::BoundingBox()
    {
    }

    BoundingBox::BoundingBox(const kvs::UnstructuredVolumeObject* object, size_t cindex)
    {
        float max = std::numeric_limits<float>::max();
        float min = - std::numeric_limits<float>::max();
     
        float xmax(min), xmin(max), ymax(min), ymin(max), zmax(min), zmin(max);

        size_t length = object->cellType();

        const kvs::UnstructuredVolumeObject::Connections& connections = object->connections();
        const kvs::UnstructuredVolumeObject::Coords& coords = object->coords();

        for ( unsigned int i = 0; i < length; i ++ )
        {
            // connections is a array of a period of Celltype
            // for Tet1 for example:
            // connections:							 0 1 2 3 4 5 6 7 3 8 10 2 ...
            //										[0 1 2 3] [4 5 6 7] [2 8 10 3] ...
            // the connection for the 3rd(cellindex==2) cell is:	    [2 8 10 3] ...
            // the first point is the 3th point in coords arry, which lietrally should be:
            // pos = cell index * Celltype (kvs cellindex starts at 0 )
            // coords array is like:				1.0 1.0 1.0 2.1 3.2 4.2 4.5 4.6 4.7 ...
            //										[1.0 1.0 1.0] [2.1 3.2 4.2] [4.5 4.6 4.7]
            // the 3rds point starts at the 6sixth position(4.5), which literally should be:
            // coords[ 3 * pos ]
            // coords[ 3 * pos + 1 ]
            // coords[ 3 * pos + 2 ]

            size_t pos = connections[ length*cindex+i ];

            // for x
            if ( coords[ 3*pos ] < xmin ) 
            {
                xmin = coords[ 3*pos ];
            }
            if( coords[ 3*pos ] > xmax )
            {
                xmax = coords[ 3*pos ];
            }

            // for y
            if ( coords[ 3*pos+1 ] < ymin ) 
            {
                ymin = coords[ 3*pos+1 ];
            }
            if( coords[ 3*pos+1 ] > ymax )
            {
                ymax = coords[ 3*pos+1 ];
            }

            // for zs
            if ( coords[ 3*pos+2 ] < zmin ) 
            {
                zmin = coords[ 3*pos+2 ];
            }
            if( coords[ 3*pos+2 ] > zmax )
            {
                zmax = coords[ 3*pos+2 ];
            }

        }
            m_bounds[0] = xmin;
            m_bounds[1] = xmax;
            m_bounds[2] = ymin;
            m_bounds[3] = ymax;
            m_bounds[4] = zmin;
            m_bounds[5] = zmax;
    }

    BoundingBox::~BoundingBox()
    {
    }

    const float* BoundingBox::bounds() const
    {
        return m_bounds;
    }
}
