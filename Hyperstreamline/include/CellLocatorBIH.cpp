#include "CellLocatorBIH.h"

namespace kvs
{
CellLocatorBIH::CellLocatorBIH()
{
    m_builder = new CellTreeBuilder();
    m_celltree = new CellTree();
    m_cache1[0] = 0;
    m_cp1 = m_cache1 + 1;
 
 }

CellLocatorBIH::CellLocatorBIH( const kvs::UnstructuredVolumeObject* volume ):
    CellLocator( volume )
{
    m_builder = new CellTreeBuilder();
    m_celltree = new CellTree();
    m_cache1[0] = 0;
    m_cp1 = m_cache1 + 1;
  
 }

CellLocatorBIH::~CellLocatorBIH()
{
    delete m_builder;
    delete m_celltree;
}

const CellTree* CellLocatorBIH::cellTree() const
{
   return m_celltree;
}

void CellLocatorBIH::setCellTree( const kvs::CellTree* ct )
{
    m_celltree = const_cast<kvs::CellTree*>(ct);
    //const_cast<const kvs::CellTree*>(m_celltree);
}

void CellLocatorBIH::setParallel()
{
	m_builder->setParallel();
}


bool CellLocatorBIH::write( const std::string filename )
{
    std::ofstream outfile( filename, std::ios_base::out|std::ios_base::binary );
    
    if (!outfile.is_open())
    {
        std::cout << "can't open " << filename << " for writing. write failed" << std::endl;
        return false;
    }

    // number of nodes
    const unsigned int nnodes = m_celltree->nodes.size();
    outfile << nnodes;

    for ( unsigned int i = 0; i < nnodes; i ++ )
    {
        outfile << m_celltree->nodes[i].index
            << m_celltree->nodes[i].lm
            << m_celltree->nodes[i].rm;
    }
   
    // number of leaves
    const unsigned int nleaves = m_celltree->leaves.size();
    outfile << nleaves;

    for ( unsigned int i = 0; i < nleaves; i ++ )
    {
        outfile << m_celltree->leaves[i];
    }
    
    outfile.close();
    return true;
}
bool CellLocatorBIH::read( const std::string filename )
{
    std::ifstream infile( filename, std::ios_base::in|std::ios_base::binary );

    if( !infile.is_open() )
    {
        std::cout << "can't open " << filename << " for reading. read failed" << std::endl;
        return false;
    }

    unsigned int nnodes;
    infile >> nnodes;
    
    kvs::CellTree* ct = new kvs::CellTree();
    ct->nodes.resize( nnodes );
    for ( unsigned int i = 0; i < nnodes; i ++ )
    {
        unsigned int index;
        float lm,rm;
        infile >> index >> lm >> rm;
        kvs::CellTree::node node;
        node.index = index;
        node.lm = lm;
        node.rm = rm;
        ct->nodes[i] = node;
    }

    unsigned int nleaves;
    infile >> nleaves;

    ct->leaves.resize( nleaves );
    for ( unsigned int i = 0; i < nleaves; i ++ )
    {
        unsigned int sz;
        infile >> sz;
        ct->leaves[i] = sz;
    }
    return true;
}
bool CellLocatorBIH::check( const kvs::CellLocator* locator )
{
    return true;
}



void CellLocatorBIH::build()
{
    if( BaseClass::m_dataset == 0 )
        std::cerr << "No dataset! \n";
    
    kvs::Timer timer;
    timer.start();
    m_builder->build( *m_celltree, BaseClass::m_dataset );
    timer.stop();
    std::cout << "Cell Tree ..." ;
    std::cout << "\t " << timer.sec() << " seconds " << std::endl;

}

int CellLocatorBIH::findCell( const float pos[3] )
{
    if ( m_celltree == 0 )
        return -1;

    switch ( BaseClass::m_mode )
    {
    case BaseClass::CACHEOFF:
    {
        CellTree::pre_traversal pt(*(this->m_celltree), pos );
        //pt.next() brings us to a series of leaves that may contain pos[3]
        while ( const CellTree::node* n = pt.next() )  
        {
            const unsigned int* begin = &(this->m_celltree->leaves[n->start()]);
            const unsigned int* end   = begin + n->size();

            for( ; begin != end; ++begin )
            {
                if ( BaseClass::testCell( *begin, pos ) )
                {
                    return *begin;
                }
            }
        }
        break;
    }
    case BaseClass::CACHEHALF:  // use cache1 and cp1
    {
        CellTree::pre_traversal_cached pt( *(this->m_celltree), pos, m_cache1[0] );
        //pt.next() brings us to a series of leaves that may contain pos[3]
        while ( const CellTree::node* n = pt.next() )  
        {
            const unsigned int* begin = &(this->m_celltree->leaves[n->start()]);
            const unsigned int* end   = begin + n->size();

            for( ; begin != end; ++begin )
            {
                if ( BaseClass::testCell( *begin, pos ) )
                {
                    const unsigned int* sp = pt.sp();
                    m_cache1[0] =  *sp;

                    return *begin;
                }
            }
        }
        break;
    }
    case BaseClass::CACHEFULL:  // use both cache stack
    {
        CellTree::in_traversal_cached pt( *(this->m_celltree), pos, m_cache1, m_cp1, m_cache2, m_cp2 );
        //pt.next() brings us to a series of leaves that may contain pos[3]
        while ( const CellTree::node* n = pt.next() )  
        {
            const unsigned int* begin = &(this->m_celltree->leaves[n->start()]);
            const unsigned int* end   = begin + n->size();

            for( ; begin != end; ++begin )
            {
                if ( BaseClass::testCell( *begin, pos ) )
                {
                    const unsigned int* stack1 = pt.stack();
                    const unsigned int* sp1 = pt.sp();
                    int n = sp1-stack1;

                    memcpy( m_cache1+1, stack1+1, 124 );
                    m_cp1 = m_cache1 + n + 1;                   //+1 is important!!

                    const unsigned int* stack2 = pt.lr_stack();
                    const unsigned int* sp2 = pt.lr_sp();
                    n = sp2-stack2;

                    memcpy( m_cache2, stack2, 64 );
                    if ( n > 1 )
                        n --;
                    m_cp2 = m_cache2 + n;                       //-1 is important!!

                    return *begin;
                }
            }
        }
        break;
    }
    default:
    {
        std::cerr << "CellLocator Mode parameter error!\n";
    }
    }


    return -1;
}
 
bool CellLocatorBIH::isDegenerate() const
{
    return false;
}

void CellLocatorBIH::clearCache()
{

    for( int i = 0; i < 32; i++ )
    {
        m_cache1[ i ] = 0;
    }

    for( int i = 0; i < 16; i++ )
    {
        m_cache2[ i ] = 0;
    }

    m_cp1 = m_cache1+1;
    m_cp2 = m_cache2;

}

}
