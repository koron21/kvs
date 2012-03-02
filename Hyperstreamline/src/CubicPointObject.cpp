#include <kvs/ValueArray>
#include "CubicPointObject.h"

namespace kvs
{

CubicPointObject::~CubicPointObject()
{
}

void CubicPointObject::reset_coordinates()
{
    std::vector<kvs::Real32> v;
    for( unsigned int i = 0; i <= m_resx; i++ )
    {
        for ( unsigned int j = 0; j <= m_resy; j ++ )
        {
            for ( unsigned int k = 0; k <= m_resz; k ++ )
            {
                v.push_back( ( m_xmax - m_xmin ) / m_resx * i + m_xmin );
                v.push_back( ( m_ymax - m_ymin ) / m_resy * j + m_ymin );
                v.push_back( ( m_zmax - m_zmin ) / m_resz * k + m_zmin );
            }
        }   
    }

    kvs::ValueArray<float> coords( v );
    kvs::ValueArray<unsigned char> colors( coords.size() );

    for ( unsigned int m = 0; m < colors.size() / 3; m ++ )
    {
        colors[ 3 * m ] = 255;
        colors[ 3 * m + 1 ] = 0;
        colors[ 3 * m + 2 ] = 0;
    }
    this->setCoords( coords );
    this->setColors( colors );
}

void CubicPointObject::reset_coordinates( 
        const unsigned int resx, const unsigned int resy, const unsigned int resz, 
        const float xmin, const float xmax, const float ymin, 
        const float ymax, const float zmin, const float zmax )
{
    this->setResX( resx );
    this->setResY( resy );
    this->setResZ( resz );
    this->setXMin( xmin );
    this->setXMax( xmax );
    this->setYMin( ymin );
    this->setYMax( ymax );
    this->setZMin( zmin );
    this->setZMax( zmax );
    this->reset_coordinates();
}

void CubicPointObject::request_update( void (*update_streamline)() )
{
    update_streamline();
}

}
