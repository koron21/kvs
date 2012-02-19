#ifndef KVS_CUBIC_POINT_OBJECT__INCLUDE
#define KVS_CUBIC_POINT_OBJECT__INCLUDE

#include <kvs/ClassName>
#include <kvs/Module>
#include <kvs/PointObject>

namespace kvs{

class CubicPointObject : public kvs::PointObject
{
    kvsClassName( kvs::CubicPointObject );

    kvsModuleCategory( Object );
    kvsModuleBaseClass( kvs::PointObject );

public:

    CubicPointObject():
        kvs::PointObject(){}
 
    CubicPointObject(
        const kvs::ValueArray<kvs::Real32>& coords,
        const kvs::ValueArray<kvs::UInt8>&  colors,
        const kvs::ValueArray<kvs::Real32>& normals,
        const kvs::ValueArray<kvs::Real32>& sizes ) :
        kvs::PointObject( coords, colors, normals, sizes ){}

    CubicPointObject(
        const kvs::ValueArray<kvs::Real32>& coords,
        const kvs::ValueArray<kvs::UInt8>&  colors,
        const kvs::ValueArray<kvs::Real32>& normals,
        const kvs::Real32                   size ) :
        kvs::PointObject( coords, colors, normals, size ){}

    CubicPointObject(
        const kvs::ValueArray<kvs::Real32>& coords,
        const kvs::RGBColor&                color,
        const kvs::ValueArray<kvs::Real32>& normals,
        const kvs::ValueArray<kvs::Real32>& sizes ) :
        kvs::PointObject( coords, color, normals, sizes ){}

    CubicPointObject(
        const kvs::ValueArray<kvs::Real32>& coords,
        const kvs::ValueArray<kvs::Real32>& normals,
        const kvs::ValueArray<kvs::Real32>& sizes ) :
        kvs::PointObject( coords, normals, sizes ){} 

    CubicPointObject(
        const kvs::ValueArray<kvs::Real32>& coords,
        const kvs::RGBColor&                color,
        const kvs::ValueArray<kvs::Real32>& normals,
        const kvs::Real32                   size ) :
        kvs::PointObject( coords, color, normals, size ){}

    CubicPointObject(
        const kvs::ValueArray<kvs::Real32>& coords,
        const kvs::ValueArray<kvs::UInt8>&  colors,
        const kvs::ValueArray<kvs::Real32>& sizes ) :
        kvs::PointObject( coords, colors, sizes ){} 

    CubicPointObject(
        const kvs::ValueArray<kvs::Real32>& coords,
        const kvs::RGBColor&                color,
        const kvs::ValueArray<kvs::Real32>& sizes ) :
        kvs::PointObject( coords, color, sizes ){}

    CubicPointObject(
        const kvs::ValueArray<kvs::Real32>& coords,
        const kvs::ValueArray<kvs::UInt8>&  colors,
        const kvs::Real32                   size ) :
        kvs::PointObject( coords, colors, size ){}

    CubicPointObject(
        const kvs::ValueArray<kvs::Real32>& coords,
        const kvs::RGBColor&                color,
        const kvs::Real32                   size ) :
        kvs::PointObject( coords, color, size ){}

    CubicPointObject(
        const kvs::ValueArray<kvs::Real32>& coords ) :
        kvs::PointObject( coords ){}

    CubicPointObject( const kvs::PointObject& other ) :
        kvs::PointObject( other ){}

    CubicPointObject( const kvs::LineObject& line ) :
        kvs::PointObject( line ){}

    CubicPointObject( const kvs::PolygonObject& polygon ) :
        kvs::PointObject( polygon ){}

    ~CubicPointObject( void );

public:

    const unsigned int resX() const { return m_resx; }
    const unsigned int resY() const { return m_resy; }
    const unsigned int resZ() const { return m_resz; }
    const float xMin() const { return m_xmin; }
    const float xMax() const { return m_xmax; }
    const float yMin() const { return m_ymin; }
    const float yMax() const { return m_ymin; }
    const float zMin() const { return m_zmax; }
    const float zMax() const { return m_zmax; }

    void setResX( const unsigned int resx ) { m_resx = resx; }
    void setResY( const unsigned int resy ) { m_resy = resy; }
    void setResZ( const unsigned int resz ) { m_resz = resz; }
    void setXMin( const float xmin ) { m_xmin = xmin; }
    void setXMax( const float xmax ) { m_xmax = xmax; }
    void setYMin( const float ymin ) { m_ymin = ymin; }
    void setYMax( const float ymax ) { m_ymax = ymax; }
    void setZMin( const float zmin ) { m_zmin = zmin; }
    void setZMax( const float zmax ) { m_zmax = zmax; }

public:

    void reset_coordinates();
    void reset_coordinates( 
            const unsigned int resx, const unsigned int resy, const unsigned int resz, 
            const float xmin, const float xmax, const float ymin, 
            const float ymax, const float zmin, const float zmax );

protected:

    unsigned int        m_resx;
    unsigned int        m_resy;
    unsigned int        m_resz;
    float               m_xmin;
    float               m_xmax;
    float               m_ymin;
    float               m_ymax;
    float               m_zmin;
    float               m_zmax;


};


}

#endif
