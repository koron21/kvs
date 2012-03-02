#ifndef KVS__CELL_LOCATOR_BIH_H_INCLUDE
#define KVS__CELL_LOCATOR_BIH_H_INCLUDE

#include "CellLocator.h"
#include "CellTree.h"

#include <iostream>
#include <fstream>
#include <kvs/ClassName>
#include <kvs/Timer>

namespace kvs
{


class CellLocatorBIH : public CellLocator
{
    kvsClassName( kvs::CellLocatorBIH )

public:

    typedef CellLocator BaseClass;

public:
        
    CellLocatorBIH();
    CellLocatorBIH( const kvs::UnstructuredVolumeObject* volume );
    ~CellLocatorBIH();

public:

    const kvs::CellTree* cellTree( void ) const;
    void setCellTree( const kvs::CellTree* ct );
    void setParallel();

public:

    void build();
    int findCell( const float pos[3] );
    bool isDegenerate() const;
    void clearCache();
    bool write( const std::string filename );
    bool read( const std::string filename );
    bool check( const kvs::CellLocator* locator );

private:

    CellTreeBuilder*                       m_builder; 
    mutable CellTree*                      m_celltree;
    unsigned int                           m_cache1[64];
    unsigned int*                          m_cp1;
    unsigned int                           m_cache2[32];
    unsigned int*                          m_cp2;

};

}

#endif
