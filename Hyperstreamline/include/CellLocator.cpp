#include "CellLocator.h"
#include <kvs/Timer>

namespace kvs
{

CellLocator::CellLocator()
{
    m_mode = CACHEOFF;
}


CellLocator::CellLocator( const kvs::UnstructuredVolumeObject* volume ):
m_dataset(volume)
{
    m_mode = CACHEOFF;
    this->initializeCell();
}

void CellLocator::initializeCell()
{
    if ( m_dataset == NULL )
        std::cerr << "no dataset!" << std::endl;
    
    switch( m_dataset->cellType() )
    {
        case kvs::UnstructuredVolumeObject::QuadraticTetrahedra:
        {
            m_cell = new kvs::QuadraticTetrahedralCell<float>(m_dataset);
            break;
        }
        case kvs::UnstructuredVolumeObject::Tetrahedra:
        {
            m_cell = new kvs::TetrahedralCell<float>(m_dataset);
            break;
        }
        default:
        {
            break;
        }
    }
}

CellLocator::~CellLocator()
{       
    delete m_cell;
}

bool CellLocator::testCell( size_t cellid, const float pos[3] ) const
{
//    static int counter = 0;
//    static float time_in_all = 0;
//    counter++;
//    kvs::Timer t;
//    t.start();

    kvs::BoundingBox bd(m_dataset, cellid);
    const float* bounds = bd.bounds();

    //if it's outside the bounding box of the cell return false
    if( pos[0]<bounds[0] || pos[0]>bounds[1] ||
        pos[1]<bounds[2] || pos[1]>bounds[3] ||
        pos[2]<bounds[4] || pos[2]>bounds[5] )
        return false;
    
    kvs::Vector3f gpos( pos[0], pos[1], pos[2] );
    m_cell->bindCell( cellid );
    m_cell->setGlobalPoint( gpos );
    kvs::Vector3f local = m_cell->localPoint();    
    //std::cout << local << std::endl;
    //

//    t.stop();
//    time_in_all += t.msec();
//    if ( counter == 72399 )
//        std::cout << "\npoint-in-cell test costs" << time_in_all << "\tmsecs" << std::endl;
    

    for ( int t = 0; t < 3; t++ )
    {
        if ( local[t] <= 0 || local[t] >= 1 )
            return false;
    }

    if ( local[0] + local[1] + local[2] > 1 )
        return false;

    return true;
}

kvs::Vector3f CellLocator::center( size_t cellid ) const
{
    m_cell->bindCell( cellid );
    int nnodes = m_dataset->cellType();

    kvs::Vector3f center;

    for ( int i = 0; i < nnodes; i ++ )
    {
        center.x() += m_cell->vertices()[i].x();
        center.y() += m_cell->vertices()[i].y();
        center.z() += m_cell->vertices()[i].z();
    }

    center.x() /= nnodes; 
    center.y() /= nnodes; 
    center.z() /= nnodes; 

    return center;
}

size_t CellLocator::randomCellindex() const
{

    return std::rand() % m_dataset->ncells();
}

void CellLocator::printVertices() const
{
    for ( int i = 0; i < m_dataset->cellType() ; i++ )
    {
        std::cout << "v[" << i << "] = " << m_cell->vertices()[i] << std::endl;
    }
}

}
