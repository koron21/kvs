#include "CellLocatorBIH.h"

namespace kvs
{
CellLocatorBIH::CellLocatorBIH()
{
    m_builder = new CellTreeBuilder();
    m_celltree = new CellTree();

    this->clearCache();
 }

CellLocatorBIH::CellLocatorBIH( const kvs::UnstructuredVolumeObject* volume ):
    CellLocator( volume )
{
    m_builder = new CellTreeBuilder();
    m_celltree = new CellTree();
  
    this->clearCache();
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
    //std::ofstream outfile( filename.c_str(), std::ios_base::out|std::ios_base::binary );
    std::ofstream outfile( filename.c_str() );
    
    if (!outfile.is_open())
    {
        std::cout << "can't open " << filename << " for writing. write failed" << std::endl;
        return false;
    }

    // number of nodes
    const unsigned int nnodes = m_celltree->nodes.size();
    outfile << nnodes << std::endl;

    for ( unsigned int i = 0; i < nnodes; i ++ )
    {
        unsigned int index = m_celltree->nodes[i].index;
        outfile << index << " ";
        if ( index ==3 )
        {
            outfile << m_celltree->nodes[i].st << " "
                    << m_celltree->nodes[i].sz << std::endl;
        }
        else
        {
            outfile << m_celltree->nodes[i].lm << " "
                    << m_celltree->nodes[i].rm << std::endl;
        }

    }
   
    // number of leaves
    const unsigned int nleaves = m_celltree->leaves.size();
    outfile << nleaves << std::endl;

    for ( unsigned int i = 0; i < nleaves; i ++ )
    {
        outfile << m_celltree->leaves[i] << std::endl;
    }
    
    outfile.close();
    return true;
}
bool CellLocatorBIH::read( const std::string filename )
{
    //std::ifstream infile( filename.c_str(), std::ios_base::in|std::ios_base::binary );
    std::ifstream infile( filename.c_str() );

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
        unsigned int st,sz;
        infile >> index;
        
        kvs::CellTree::node node;
        node.index = index;

        if ( index == 3 )
        {
            infile >> st >> sz;
            node.st = st;
            node.sz = sz;
        }
        else
        {
            infile >> lm >> rm;
            node.lm = lm;
            node.rm = rm;
        }

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

    if ( m_celltree!= NULL )
    {
        delete m_celltree;
        m_celltree = NULL;
    }
    m_celltree = ct;
    return true;
}
bool CellLocatorBIH::check( const kvs::CellLocator* locator )
{
    const kvs::CellLocatorBIH* bih = static_cast<const kvs::CellLocatorBIH*>( locator );
    const kvs::CellTree* ct1 = bih->cellTree();
    const kvs::CellTree* ct2 = m_celltree;
    
    //check nodes
    const unsigned int node_size1 = ct1->nodes.size();
    const unsigned int node_size2 = ct2->nodes.size();

    if ( node_size1 != node_size2 )
        return false;

    for ( unsigned int i = 0; i < node_size1; i ++ )
    {
        bool flag1 = ( ct1->nodes[i].index == ct2->nodes[i].index );
        if ( !flag1 )
        {
            std::cout << "\ntree1 and tree2 differs at " << i << "th node's index. " << std::endl;
            return false;
        }
        
        bool flag2 = false;
        if ( ct1->nodes[i].index == 3 ) // if it is a leaf
        {
            flag2 = ( ct1->nodes[i].st == ct2->nodes[i].st    ) &&
                    ( ct1->nodes[i].sz == ct2->nodes[i].sz    );
        }
        else
        {
            flag2 = ( std::abs( ct1->nodes[i].lm - ct2->nodes[i].lm ) < 0.001    ) &&
                    ( std::abs( ct1->nodes[i].rm - ct2->nodes[i].rm ) < 0.001   );
        }

        if ( !flag2 )
        {
            std::cout << "\ntree1 and tree2 differs at " << i << "th node. " << std::endl;
            std::cout << "Tree1, node" << i << std::endl;
            std::cout << "Index      :" << ct1->nodes[i].index << std::endl;
            if ( ct1->nodes[i].index != 3 )
            {
                std::cout << "lmax       :" << ct1->nodes[i].lm << std::endl;
                std::cout << "Rmin       :" << ct1->nodes[i].rm << std::endl;
            }
            else
            {
                std::cout << "Start      :" << ct1->nodes[i].st << std::endl;
                std::cout << "Size       :" << ct1->nodes[i].sz << std::endl;
            }
            std::cout << "Tree2, node" << i << std::endl;
            std::cout << "Index      :" << ct2->nodes[i].index << std::endl;
            if ( ct2->nodes[i].index != 3 )
            {
                std::cout << "lmax       :" << ct2->nodes[i].lm << std::endl;
                std::cout << "Rmin       :" << ct2->nodes[i].rm << std::endl;
            }
            else
            {
                std::cout << "Start      :" << ct2->nodes[i].st << std::endl;
                std::cout << "Size       :" << ct2->nodes[i].sz << std::endl;
            }
            return false;
        }
    }
    
    //check leaves
    const unsigned int leaf_size1 = ct1->leaves.size();
    const unsigned int leaf_size2 = ct2->leaves.size();

    if ( leaf_size1 != leaf_size2 )
        return false;

    for ( unsigned int i = 0; i < leaf_size1; i ++ )
    {
        bool flag = ( ct1->leaves[i] == ct2->leaves[i] );

        if (!flag)
        {
            std::cout << "tree1 and tree2 differs at " << i << "th leaf. " << std::endl;
            return false;
        }
    }

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
        CellTree::in_traversal_cached pt( *(this->m_celltree), pos, m_cache1, m_cp1 );
        //pt.next() brings us to a series of leaves that may contain pos[3]
        while ( const CellTree::node* n = pt.next() )  
        {
            const unsigned int* begin = &(this->m_celltree->leaves[n->start()]);
            const unsigned int* end   = begin + n->size();

            for( ; begin != end; ++begin )
            {
                if ( BaseClass::testCell( *begin, pos ) )
                {
                    const unsigned int* stack1 = pt.m_stack;
                    const unsigned int* sp1 = pt.m_sp;
                    int n = sp1-stack1;

                    memcpy( m_cache1+1, stack1+1, 124 );
                    m_cp1 = m_cache1 + n + 1;                   //+1 is important!!

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
        m_cache1[ i ] = -1;
    }

    for( int i = 0; i < 16; i++ )
    {
        m_cache2[ i ] = -1;
    }

    m_cache1[0] = 0;
    m_cache2[0] = 0;

    m_cp1 = m_cache1+1;
    m_cp2 = m_cache2;

}

}
