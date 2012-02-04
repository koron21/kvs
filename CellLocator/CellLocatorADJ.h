#ifndef KVS__CELL_LOCATOR_ADJ_H_INCLUDE
#define KVS__CELL_LOCATOR_ADJ_H_INCLUDE

#include "CellLocator.h"
#include <kvs/CellAdjacencyGraph>
#include <kvs/ClassName>
#include <kvs/Timer>

namespace kvs
{

class CellLocatorADJ : public CellLocator
{
    kvsClassName(kvs::CellLocatorADJ);

public:

    typedef kvs::CellLocator BaseClass;

public:

    CellLocatorADJ();
    CellLocatorADJ( const kvs::UnstructuredVolumeObject* volume );
    ~CellLocatorADJ();

public:
    const kvs::CellAdjacencyGraph* adjacencyGraph() const;

public:

    void build();
    int findCell( const float pos[3] );
    int find( const float pos[3], const int start_cellid );
    bool isDegenerate() const;
    void clearCache();

    // write the datasturcture in to binary file
    virtual bool write( const std::string filename );
    // read the datastructure from binary file
    virtual bool read( const std::string filename );
    // check whether two locator generates the same datastructure
    virtual bool check( const kvs::CellLocator* locator );

private:

    const kvs::CellAdjacencyGraph*  m_adj;
    unsigned int                    m_nrandtests;
    int                             m_hint_cellid;   // used for cache
};


}
#endif
