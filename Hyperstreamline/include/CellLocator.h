#ifndef KVS__CELL_LOCATOR_H_INCLUDE
#define KVS__CELL_LOCATOR_H_INCLUDE

#include "BoundingBox.h"

#include <kvs/ClassName>
#include <kvs/CellBase>
#include <kvs/UnstructuredVolumeObject>
#include <kvs/QuadraticTetrahedralCell>
#include <kvs/TetrahedralCell>
#include <kvs/Vector>
#include <kvs/Matrix>

namespace kvs{


class CellLocator
{
    kvsClassName( kvs::CellLocator );

public:
    
    enum Mode
    {
        CACHEOFF = 0,
        CACHEHALF  = 1,
        CACHEFULL = 2
    };

public:

    CellLocator();
    CellLocator( const kvs::UnstructuredVolumeObject* volume );
    virtual ~CellLocator();

public:

    bool testCell( size_t cellid, const float pos[3] ) const;
    kvs::Vector3f center( size_t cellid ) const;
    size_t randomCellindex() const;
    void printVertices() const;

    virtual void build() = 0;
    virtual int findCell( const float pos[3] ) = 0;
    virtual bool isDegenerate() const = 0;
    virtual void clearCache() = 0;
    // write the datasturcture in to binary file
    virtual bool write( const std::string filename ) = 0;
    // read the datastructure from binary file
    virtual bool read( const std::string filename ) = 0;
    // check whether two locator generates the same datastructure
    virtual bool check( const kvs::CellLocator* locator ) = 0;

public:

    const CellLocator::Mode mode() const { return m_mode; } 
    const kvs::UnstructuredVolumeObject* dataset() { return m_dataset; }
    const kvs::CellBase<float>* cell() { return m_cell; }
    void setDataSet( const kvs::UnstructuredVolumeObject* volume ){ m_dataset = volume; }
    void initializeCell();
    void setMode( const CellLocator::Mode mode ){ m_mode = mode; }

protected:

    const kvs::UnstructuredVolumeObject*    m_dataset;
    kvs::CellBase<float>*                   m_cell;
    Mode                                    m_mode;
    
};

}
#endif
