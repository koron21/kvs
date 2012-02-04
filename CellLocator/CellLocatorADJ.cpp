#include "CellLocatorADJ.h"
#include <kvs/Matrix>


namespace {

struct Line
{
    Line() {}
    Line( const kvs::Vector3f st, const kvs::Vector3f ed, const float T ) :
        start(st), end(ed), t(T)
    {
    }
    kvs::Vector3f start;
    kvs::Vector3f end;
    float         t;
};

struct Plane
{
    Plane() {}
    Plane( const kvs::Vector3f a, const kvs::Vector3f b, const kvs::Vector3f c ) :
        p0(a), p1(b), p2(c)
    {
    }

    kvs::Vector3f p0;
    kvs::Vector3f p1;
    kvs::Vector3f p2;
};

struct weights  //the return value of linePlaneIntersection
{
    float t;
    float u;
    float v;
    bool parallel;
};

const kvs::UInt32 TetrahedralCellFaces[12] = 
{
    0, 1, 2, // face 0
    0, 2, 3, // face 1
    0, 3, 1, // face 2
    1, 3, 2  // face 3
};

const float tiny_value = 0.1;

}

namespace kvs
{

CellLocatorADJ::CellLocatorADJ()
{
    m_hint_cellid = -1;
    m_nrandtests = 30;
}

CellLocatorADJ::CellLocatorADJ( const kvs::UnstructuredVolumeObject* volume ) : CellLocator( volume )
{
    m_hint_cellid = -1;
    m_nrandtests = 30;
}

CellLocatorADJ::~CellLocatorADJ()
{
    delete m_adj;
}

const kvs::CellAdjacencyGraph* CellLocatorADJ::adjacencyGraph() const
{
    return m_adj;
}

void CellLocatorADJ::build()
{
    if ( !BaseClass::m_dataset )
    {
        std::cerr << "no data to build adjacency table!!\n" << std::endl;
        exit(-1);
    }
    kvs::Timer timer;
    timer.start();
    std::cout << "Building Adjacency Graph ... ";
    m_adj = new kvs::CellAdjacencyGraph( BaseClass::m_dataset );
    timer.stop();
    std::cout << timer.sec() << " seconds" << std::endl;
}

bool CellLocatorADJ::isDegenerate() const
{
    switch ( BaseClass::m_dataset->cellType() )
    {
        case kvs::UnstructuredVolumeObject::Tetrahedra:
        {
            kvs::Vector3f V1 = BaseClass::m_cell->vertices()[0] - BaseClass::m_cell->vertices()[0];
            kvs::Vector3f V2 = BaseClass::m_cell->vertices()[1] - BaseClass::m_cell->vertices()[0];
            kvs::Vector3f V3 = BaseClass::m_cell->vertices()[2] - BaseClass::m_cell->vertices()[0];

            kvs::Matrix33f M( V1, V2, V3 );
            if ( M.determinant() == 0 )
                return true;
            else
                return false;
        }
        default:
        {
            return false;
        }
    }
}

const weights linePlaneIntersection( Line line, const Plane plane )
{
    weights w;
    const kvs::Vector3f& la = line.start;
    const kvs::Vector3f& lb = line.end;
    //const float& t          = line.t;
    
    const kvs::Vector3f& p0 = plane.p0; 
    const kvs::Vector3f& p1 = plane.p1; 
    const kvs::Vector3f& p2 = plane.p2; 

    kvs::Matrix33f M(
            la.x()-lb.x(),  p1.x()-p0.x(),  p2.x()-p0.x(),
            la.y()-lb.y(),  p1.y()-p0.y(),  p2.y()-p0.y(),
            la.z()-lb.z(),  p1.z()-p0.z(),  p2.z()-p0.z()
            );

    if ( M.determinant()==0 ) // if singular
    {
        w.t = 0;
        w.u = 0;
        w.v = 0;
        w.parallel = 1;
        return w;
    }

    kvs::Vector3f V(
            la.x()-p0.x(),
            la.y()-p0.y(),
            la.z()-p0.z()
            );

    kvs::Vector3f result = M.inverse() * V;
    w.t = result[0];
    w.u = result[1];
    w.v = result[2];
    w.parallel = 0;

    return w;
}

int CellLocatorADJ::findCell( const float pos[3] ) 
{
    switch ( BaseClass::m_mode )
    {
        case CellLocator::CACHEOFF:
        {
            unsigned int startindex, temp_startindex;
            float min = std::numeric_limits<float>::max();
            float distance;
            kvs::Vector3f center;

            // 1 bind some random cell indices, get their center, find the one closest to the target point
            for ( unsigned int i = 0; i < m_nrandtests; i ++ )
            {
                temp_startindex = this->randomCellindex();
                BaseClass::m_cell->bindCell( temp_startindex );
                center = this->center( temp_startindex );
                distance = std::pow( center.x()-pos[0], 2 ) + std::pow( center.y()-pos[1], 2 ) + std::pow( center.z()-pos[2], 2 );
                if ( distance < min )
                {
                    min = distance;
                    startindex = temp_startindex;
                }
            }

            return ( this->find( pos, startindex ) );
        }
        case CellLocator::CACHEHALF:
        {
            if ( m_hint_cellid == -1 )
            {   
                // 1 bind some random cell indices, get their center, find the one closest to the target point
                unsigned int startindex, temp_startindex;
                float min = std::numeric_limits<float>::max();
                float distance;
                kvs::Vector3f center;

                for ( unsigned int i = 0; i < m_nrandtests; i ++ )
                {
                    temp_startindex = this->randomCellindex();
                    BaseClass::m_cell->bindCell( temp_startindex );
                    center = this->center( temp_startindex );
                    distance = std::pow( center.x()-pos[0], 2 ) + std::pow( center.y()-pos[1], 2 ) + std::pow( center.z()-pos[2], 2 );
                    if ( distance < min )
                    {
                        min = distance;
                        startindex = temp_startindex;
                    }
                }

                m_hint_cellid = startindex;     
            }

            return ( this->find( pos, m_hint_cellid ) );
        }
        default:
        {
            std::cerr << "CellLocator::Mode receives wrong parameter!\n" ;
        }
    }

}

int CellLocatorADJ::find( const float pos[3], const int start_cellid )
{
    switch ( BaseClass::m_dataset->cellType())
    {
        case kvs::UnstructuredVolumeObject::QuadraticTetrahedra:
        case kvs::UnstructuredVolumeObject::Tetrahedra:
        {
            unsigned int current_faceid; 
            weights w;
            float step = 0;
            bool found = false;

            // 1 starting from the center, find the intersection of the line and polygon
            //
            // 2 from adjacency graph, find which cell to go next
            // 3 go to the next cell, find the outgoing intersection
            //
            // repeat from 2 to 3 util reach the pos

            BaseClass::m_cell->bindCell( start_cellid );
            kvs::Vector3f center = this->center( start_cellid );
            kvs::Vector3f end( pos[0], pos[1], pos[2] );

            //std::cout.width( 30 );
            //std::cout << std::left << "start from" << center << std::endl;
            Line line( center, end, 0 ); //initialize the line

            int current_cellid = start_cellid;

            while (1)
            {
                if (found)
                    return current_cellid;

                //std::cout << "step = " << step << std::endl;
                if ( BaseClass::testCell( current_cellid, pos ) )
                {
                    found = true;
                    m_hint_cellid = current_cellid;
                    return current_cellid;
                }

                BaseClass::m_cell->bindCell( current_cellid );
                for( int i = 0; i < 4; i ++ )
                {
                    Plane p( BaseClass::m_cell->vertices()[TetrahedralCellFaces[ 3*i    ]], 
                             BaseClass::m_cell->vertices()[TetrahedralCellFaces[ 3*i+1  ]],
                             BaseClass::m_cell->vertices()[TetrahedralCellFaces[ 3*i+2  ]]);

                    w = linePlaneIntersection( line, p );
                    if ( w.u >= 0 && w.v >= 0 && w.u + w.v <= 1 && w.t > step )
                    {
                        current_faceid = i;
                        step = w.t;
                        break;
                    }

                    if ( i == 3 )
                    {
                        step = 0;
                        line.start = BaseClass::m_cell->randomSampling();
                        i = 0;
                        //step = 0;
                        //current_cellid = this->findCell( pos );
                        //this->printVertices();
                        //std::cout << "Degenerate?? : " << this->isDegenerate() << std::endl;
                        //exit(0);
                    }
                }
                if ( m_adj->mask()[ 4*current_cellid+current_faceid ] == 1 ) 
                {
                    current_cellid = m_adj->graph()[ 4*current_cellid+current_faceid ];
                }
                else // the ray goes out of the volume object
                {
                    // record current step
                    // do a brute force search along all the external face
                    // find out the step that satisfies
                    // step < 1 and step = max of all the step found
                    // and then set the current index 
                    float step_save = step;

                    for ( unsigned int i = 0; i < this->m_adj->mask().size(); i ++ )
                    {
                        if ( m_adj->mask()[i] == 0 )
                        {
                            current_faceid = i % 4;
                            BaseClass::m_cell->bindCell( i / 4 ); //current_cellid = i / 4;
                            Plane p( BaseClass::m_cell->vertices()[ TetrahedralCellFaces[ 3*current_faceid ] ],
                                BaseClass::m_cell->vertices()[ TetrahedralCellFaces[ 3*current_faceid+1 ] ],
                                BaseClass::m_cell->vertices()[ TetrahedralCellFaces[ 3*current_faceid+2 ] ]);
                            w = linePlaneIntersection( line, p );

                            if ( w.u >= 0 && w.v >= 0 && w.u + w.v <= 1 && w.t > step && w.t < 1 )
                            {
                                current_cellid = i / 4;
                                step = w.t;
                            }
                        }
                    }

                    if ( step_save == step ) //means the point is outside of the volume
                        return -1;
                }
            }
        }
        default:
        {
            return -1;
        }
    }
}

void CellLocatorADJ::clearCache()
{
    m_hint_cellid = -1;
}

// write the datasturcture in to binary file
bool CellLocatorADJ::write( const std::string filename )
{
    return false;
}

// read the datastructure from binary file
bool CellLocatorADJ::read( const std::string filename )
{
    return false;
}

// check whether two locator generates the same datastructure
bool CellLocatorADJ::check( const kvs::CellLocator* locator )
{
    return false;
}

}

