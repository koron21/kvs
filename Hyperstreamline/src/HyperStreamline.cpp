#include <kvs/BitArray>
#include <kvs/EigenDecomposer>
#include <kvs/IgnoreUnusedVariable>
#include <kvs/Math>
#include <kvs/Message>
#include <kvs/RGBColor>
#include <kvs/RGBFormulae>
#include <kvs/SVSolver>
#include <kvs/Thread>
#include <kvs/Type>
#include <kvs/Vector3>
#include <kvs/VolumeObjectBase>

#include <cmath>
#include "HyperStreamline.h"

#define EIGEN_VALUE_COLOR
//#define MISES_COLOR

namespace kvs
{

HyperStreamline::HyperStreamline( void ):
    kvs::MapperBase(),
    m_seed_points( NULL ),
    m_integration_method( HyperStreamline::RungeKutta2nd ),
    m_integration_direction( HyperStreamline::BothDirections ),
    m_integration_interval( 0.35f ),
    m_vector_length_threshold( 0.000001f ),
    m_integration_times_threshold( 1024 ),
    m_enable_boundary_condition( true ),
    m_enable_vector_length_condition( false ),
    m_enable_integration_times_condition( true ),
    m_degenerate( false ),
    m_nthreads( 0 ),
    m_enable_cache( true )
{
}

HyperStreamline::HyperStreamline(
    const kvs::VolumeObjectBase* volume,
    const kvs::PointObject* seed_points,
    const kvs::TransferFunction& transfer_function ):
    m_integration_interval( 0.35f ),
    m_vector_length_threshold( 0.000001f ),
    m_integration_times_threshold( 1024 ),
    m_enable_boundary_condition( true ),
    m_enable_vector_length_condition( false ),
    m_enable_integration_times_condition( true ),
    m_degenerate( false ),
    m_nthreads( 0 ),
    m_enable_cache( true )
{
    BaseClass::setTransferFunction( transfer_function );
    setSeedPoints( seed_points );
    this->exec( volume );
}

HyperStreamline::~HyperStreamline( void )
{
    m_eigenvalues.clear();
}

HyperStreamline::SuperClass* HyperStreamline::exec( const kvs::ObjectBase* object )
{
    if ( !object )
    {
        BaseClass::m_is_success = false;
        kvsMessageError("Input object is NULL.");
        return( NULL );
    }

    const kvs::VolumeObjectBase* volume = kvs::VolumeObjectBase::DownCast( object );
    if ( !volume )
    {
        BaseClass::m_is_success = false;
        kvsMessageError("Input object is not volume dat.");
        return( NULL );
    }

    if ( !volume->hasMinMaxValues() )
        volume->updateMinMaxValues();
    
    BaseClass::attach_volume( volume );
    mapping(volume);
    
    return (this);
}

void HyperStreamline::mapping( const kvs::VolumeObjectBase* volume )
{
    
    kvs::Timer timer;
    timer.start();
    std::cout << "Extracting Hyper Streamline";
    this->extract_lines( volume );
    timer.stop();
    std::cout << "\t" << std::setw( 10 ) << std::left << std::setprecision( 5 ) 
        << timer.sec() << " seconds\t" << "nVertices:\t" << this->nvertices() << std::endl;
}

void HyperStreamline::extract_lines(
    const kvs::VolumeObjectBase* volume )
{
    kvs::IgnoreUnusedVariable( volume );

    // Calculated data arrays.
    std::vector<kvs::Real32> coords;
    std::vector<kvs::UInt32> connections;
    std::vector<kvs::UInt8>  colors;

    //for each seed point, initiates a thread to calculate streamline
    const size_t npoints = m_seed_points->nvertices();
    kvs::BitArray mask( npoints );
    mask.reset();
    for ( size_t index = 0; index < npoints; index++ )
    {
        const kvs::Vector3f seed_point = m_seed_points->coord( index );
        if ( this->check_for_inside_volume( seed_point ) )
        {
            m_nthreads ++;
            mask.set(index);
        }
    }

#ifdef DEBUG
    const size_t nseeds = m_seed_points->nvertices();
    const float* p_seeds_coords = m_seed_points->coords().pointer();

    std::cout << std::endl;
    std::cout << "number of seed points:   " << nseeds << std::endl;
    for ( size_t i = 0; i < nseeds; i ++ )
    {
        printf( "%-5.2f  %-5.2f  %-5.2f  %-4lu\n", *( p_seeds_coords + 3 * i ), 
            *( p_seeds_coords + 3 * i + 1 ), *( p_seeds_coords + 3 * i + 2 ), i );
    }

    std::cout << std::endl
              << "number of worker thread: " << m_nthreads << std::endl;
    std::cout << "mask";
    for ( size_t i = 0, j = 0; i < mask.size(); i ++, j++ )
    {
        if ( j > 1 && j % 10 == 0 )
            std::cout << " ";
        if ( j % 50 == 0 ) 
            std::cout << std::endl;
        std::cout << mask[i];
    }
    std::cout << std::endl << std::endl;
#endif

    if ( m_nthreads == 0 )
    {
        coords.push_back( 0.0 );
        coords.push_back( 0.0 );
        coords.push_back( 0.0 );
        connections.push_back( 0 );
        connections.push_back( 0 );
        colors.push_back( 0 );
        colors.push_back( 0 );
        colors.push_back( 0 );
        m_eigenvalues.push_back( 1 );

        SuperClass::setLineType( kvs::LineObject::Polyline );
        SuperClass::setColorType( kvs::LineObject::VertexColor );
        SuperClass::setCoords( kvs::ValueArray<kvs::Real32>( coords ) );
        SuperClass::setConnections( kvs::ValueArray<kvs::UInt32>( connections ) );
        SuperClass::setColors( kvs::ValueArray<kvs::UInt8>( colors ) );
        SuperClass::setSize( 1.0f );
        return;
    }
    m_p_threads = new HyperStreamlineThread[m_nthreads];

    // need locator or not
    switch ( volume->volumeType() )
    {
    case kvs::VolumeObjectBase::Unstructured:
        {
            // initialize locator and set seed point
            const kvs::UnstructuredVolumeObject* uvolume = 
                static_cast<const kvs::UnstructuredVolumeObject* >(volume);

            if ( !m_locator_initialized )
            {
                m_locator = new CellLocatorBIH();
                m_locator->setDataSet( uvolume );
                m_locator->initializeCell();
                m_locator->setMode( kvs::CellLocator::CACHEOFF );
                m_locator->build();
            }
            
            for ( size_t i = 0, j = 0; (i < m_nthreads && j < npoints) ; i ++, j++ )
            {
                while ( !mask[j] )
                    j++;
                if ( j > npoints )
                    break;
                m_p_threads[i].setSeedPoint( m_seed_points->coord( j ) );
                m_p_threads[i].setLocator( m_locator->cellTree(), uvolume );    
                if ( !m_enable_cache )
                    m_p_threads[i].setDisableCache();
#ifdef DEBUG
                std::cout << j << "th seed point passed to " << i << "th thread" << std::endl;
#endif
            }

            break;
        }
    case kvs::VolumeObjectBase::Structured:
        {
            break;
        }
    default:
        {
            break;
        }
    }

#ifdef DEBUG
    std::cout << std::endl;
    std::cout << "m_nth_egvector:           " << m_nth_egvector << std::endl;
    std::cout << "m_integration_direction:  " << m_integration_direction << std::endl;
    std::cout << "m_integration_method:     " << m_integration_method << std::endl;
#endif

    //pass parameters to HyperStreamlineThread class
       for ( size_t i = 0; i < m_nthreads; i ++ )
    {
        m_p_threads[i].setGoWithNthEigenVector( m_nth_egvector );
        m_p_threads[i].setTransferFunction( BaseClass::transferFunction() );
        m_p_threads[i].setVolume( volume );
        m_p_threads[i].setDirection( this->m_integration_direction );
        m_p_threads[i].setMethod( this->m_integration_method );
        m_p_threads[i].init();
    }

    for ( size_t i = 0; i < m_nthreads; i ++ )
        m_p_threads[i].start();

    for ( size_t i = 0; i < m_nthreads; i ++ )
        m_p_threads[i].wait();

    // gather data from every thread
    size_t coords_size = 0;

    for ( size_t i = 0; i < m_nthreads; i ++ )
    {
        coords_size += m_p_threads[i].line().nvertices();
    }

#ifdef DEBUG
    std::cout << "\nnumber of vertice       " << coords_size << std::endl;
#endif

    coords.reserve( coords_size * 3 );
    colors.reserve( coords_size * 3 );
    connections.clear();
    connections.reserve( m_nthreads * 2 );
        
    for ( size_t i = 0; i < m_nthreads; i ++ )
    {
        const kvs::Real32* p_coords = m_p_threads[i].line().coords().pointer();
#ifdef MISES_COLOR
        const kvs::UInt8*  p_colors = m_p_threads[i].line().colors().pointer();
#endif

        connections.push_back( coords.size() / 3 );         //start

        for ( size_t j = 0; j < m_p_threads[i].line().nvertices(); j ++ ) 
        {
            coords.push_back( *( p_coords + 3*j ) );        //x
            coords.push_back( *( p_coords + 3*j + 1 ) );    //y
            coords.push_back( *( p_coords + 3*j + 2 ) );    //z
        
#ifdef MISES_COLOR
            colors.push_back( *( p_colors + 3*j ) );        //R
            colors.push_back( *( p_colors + 3*j + 1) );     //G
            colors.push_back( *( p_colors + 3*j + 2) );     //B
#endif
        }
        connections.push_back( coords.size() / 3 - 1 );     //end
    }
    
#ifdef EIGEN_VALUE_COLOR
    for ( size_t i = 0; i < m_nthreads; i ++ )
    {
        m_eigenvalues.insert( m_eigenvalues.end(), m_p_threads[i].eigenValue().begin(), m_p_threads[i].eigenValue().end() );
    }
    //this->calculate_color( min_eigenvalue, max_eigenvalue );
#endif

    SuperClass::setLineType( kvs::LineObject::Polyline );
    SuperClass::setColorType( kvs::LineObject::VertexColor );
    SuperClass::setCoords( kvs::ValueArray<kvs::Real32>( coords ) );
    SuperClass::setConnections( kvs::ValueArray<kvs::UInt32>( connections ) );
#ifndef EIGEN_VALUE_COLOR
    SuperClass::setColors( kvs::ValueArray<kvs::UInt8>( colors ) );
#endif
    SuperClass::setSize( 1.0f );

    delete [] m_p_threads;

}

void HyperStreamline::calculate_color()
{
    float max_eigenvalue = -std::numeric_limits<float>::max();
    float min_eigenvalue = std::numeric_limits<float>::max();

    max_eigenvalue = *std::max_element( m_eigenvalues.begin(), m_eigenvalues.end() );
    min_eigenvalue = *std::min_element( m_eigenvalues.begin(), m_eigenvalues.end() );

    this->calculate_color( min_eigenvalue, max_eigenvalue );
}

void HyperStreamline::calculate_color( const float min_eigenvalue, const float max_eigenvalue )
{
    kvs::ColorMap cmap;
    //cmap.addPoint( 0.0f, kvs::RGBColor( 0 , 155 , 0 ) ); // green
    //if ( max_eigenvalue < 0 )
    //{
    //    cmap.setRange( min_eigenvalue, 0 );
    //    cmap.addPoint( min_eigenvalue, kvs::RGBColor( 0, 0, 155 ) ); // blue 
    //    cmap.addPoint( min_eigenvalue * 7 / 8, kvs::RGBColor( 0, 0, 255 ) );  
    //    cmap.addPoint( min_eigenvalue / 8, kvs::RGBColor( 0, 255, 255 ) );  
    //}
    //if ( min_eigenvalue > 0 )
    //{
    //    cmap.setRange( 0, max_eigenvalue );
    //    cmap.addPoint( max_eigenvalue / 8, kvs::RGBColor( 55 , 155 , 0 ) ); 
    //    cmap.addPoint( max_eigenvalue * 6 / 8, kvs::RGBColor( 255 , 255 , 0 ) ); // yellow
    //    cmap.addPoint( max_eigenvalue, kvs::RGBColor( 155, 0, 0 ) ); // red
    //}
    //else
    //{
    //    cmap.setRange( min_eigenvalue, max_eigenvalue );
    //    cmap.addPoint( min_eigenvalue, kvs::RGBColor( 0 , 0 , 155 ) ); // blue
    //    cmap.addPoint( min_eigenvalue * 7 / 8, kvs::RGBColor( 0, 0, 255 ) );  
    //    cmap.addPoint( min_eigenvalue / 8, kvs::RGBColor( 0 , 255 , 255 ) ); 
    //    cmap.addPoint( max_eigenvalue / 8, kvs::RGBColor( 55 , 155 , 0 ) ); 
    //    cmap.addPoint( max_eigenvalue * 6 / 8, kvs::RGBColor( 255 , 255 , 50 ) ); // yellow
    //    cmap.addPoint( max_eigenvalue, kvs::RGBColor( 155, 0, 0 ) ); // red
    //}
    cmap.create();

    float interval;

    if ( max_eigenvalue < 0 )
    {
        interval = min_eigenvalue;
    }
    else if ( min_eigenvalue < 0 )
    {
        if ( max_eigenvalue + min_eigenvalue > 0 )
            interval = max_eigenvalue;
        else
            interval = min_eigenvalue;
    }
    else
    {
        interval = max_eigenvalue;
    }

    std::vector<kvs::UInt8>  colors;
    colors.reserve( m_eigenvalues.size()*3 );

    for ( size_t i = 0; i < m_eigenvalues.size(); i ++ )
    {
        kvs::UInt8 level;
        if ( m_eigenvalues[i] > 0 )
        {
            level = 127.0 + 128.0 * m_eigenvalues[i] / interval; // level 127.0 is zero
        }
        else
        {
            level = 128.0 * m_eigenvalues[i] / interval; // m_eigenvalues[i] < 0
        }
        colors.push_back( cmap[level].r() );
        colors.push_back( cmap[level].g() );
        colors.push_back( cmap[level].b() );
    }

    SuperClass::setColors( kvs::ValueArray<kvs::UInt8>( colors ) );
}

const bool HyperStreamline::check_for_inside_volume( const kvs::Vector3f& point )
{
    switch( BaseClass::volume()->volumeType() )
    {
        case kvs::VolumeObjectBase::Structured:
        {
            const kvs::StructuredVolumeObject* structured_volume =
            reinterpret_cast<const kvs::StructuredVolumeObject*>( BaseClass::volume() );

            const float dimx = static_cast<float>( structured_volume->resolution().x() - 1 );
            const float dimy = static_cast<float>( structured_volume->resolution().y() - 1 );
            const float dimz = static_cast<float>( structured_volume->resolution().z() - 1 );

            if ( point.x() < 0.0f || dimx < point.x() ) return( false );
            if ( point.y() < 0.0f || dimy < point.y() ) return( false );
            if ( point.z() < 0.0f || dimz < point.z() ) return( false );

            return( true );
        }
        case kvs::VolumeObjectBase::Unstructured:
        {
            if ( point.x() < this->volume()->minExternalCoord().x() || point.x() > this->volume()->maxExternalCoord().x() ) return false;
            if ( point.y() < this->volume()->minExternalCoord().y() || point.y() > this->volume()->maxExternalCoord().y() ) return false;
            if ( point.z() < this->volume()->minExternalCoord().z() || point.z() > this->volume()->maxExternalCoord().z() ) return false;

            float pos[3] = { point.x(), point.y(), point.z() };
            if ( m_locator->findCell( pos ) != -1 )
                return true;
            else
                return false;
        }
    }

    return false;

}

const std::vector<float>& HyperStreamline::eigenValues() const
{
    return m_eigenvalues;
}

void HyperStreamline::setEigenValues( const std::vector<float>& eigen_values )
{
    m_eigenvalues = eigen_values;
}

void HyperStreamline::setSeedPoints( const kvs::PointObject* seed_points )
{
    // Shallow copy.
    m_seed_points = new kvs::PointObject( seed_points->coords() );
    if ( !m_seed_points )
    {
        kvsMessageError( "Cannot allocate memory for the seed points." );
    }
}

void HyperStreamline::setIntegrationMethod( const HyperStreamline::IntegrationMethod method )
{
    m_integration_method = method;
}

void HyperStreamline::setIntegrationDirection( const HyperStreamline::IntegrationDirection direction )
{
    m_integration_direction = direction;
}

void HyperStreamline::setIntegrationInterval( const float interval )
{
    m_integration_interval = interval;
}

void HyperStreamline::setVectorLengthThreshold( const float length )
{
    m_vector_length_threshold = length;
}

void HyperStreamline::setIntegrationTimesThreshold( const size_t times )
{
    m_integration_times_threshold = times;
}

void HyperStreamline::setEnableBoundaryCondition( const bool enabled )
{
    m_enable_boundary_condition = enabled;
}

void HyperStreamline::setEnableVectorLengthCondition( const bool enabled )
{
    m_enable_vector_length_condition = enabled;
}

void HyperStreamline::setEnableIntegrationTimesCondition( const bool enabled )
{
    m_enable_integration_times_condition = enabled;
}

void HyperStreamline::setLocator( const kvs::CellTree* ct, const kvs::UnstructuredVolumeObject* volume )
{

    m_locator = new kvs::CellLocatorBIH();
    m_locator->setCellTree( ct );
    m_locator->setDataSet( volume );
    m_locator->setMode( kvs::CellLocator::CACHEOFF );
    m_locator->initializeCell();
    m_locator_initialized = true;

}

void HyperStreamline::setLocatorToInitialized()
{
    m_locator_initialized = true;
}

void HyperStreamline::setLocationMethod( const HyperStreamline::LocationMethod method, kvs::CellTree* ct, kvs::UnstructuredVolumeObject* volume )
{
    m_location_method = method;
    this->setLocator( ct, volume );
}

void HyperStreamline::setGoWithNthEigenVector( const int nth )
{
    if ( nth < 0 || nth > 2 )
    {
        std::cout << "From HyperStreamline::setGoWithNthEigenVector: There are only 3 eigen vectors!! set to default(0) " << std::endl;
        m_nth_egvector = 0;
        return;
    }
    m_nth_egvector = nth;
}

void HyperStreamline::setDisableCache()
{
    m_enable_cache = false;
}



HyperStreamlineThread::HyperStreamlineThread() : kvs::Thread()
{
    m_nth_egvector = 0;
    m_eigenvalues.reserve( 500 );
}

HyperStreamlineThread::~HyperStreamlineThread()
{
    m_eigenvalues.clear();
}

void HyperStreamlineThread::init()
{

    m_integration_interval = 0.35f;
    m_vector_length_threshold = 0.0001f;
    m_integration_times_threshold = 1024;
    m_enable_boundary_condition = true;
    m_enable_vector_length_condition = false;
    m_enable_integration_times_condition = true;
    m_degenerate = false;

}

const bool HyperStreamlineThread::check()
{
    if ( m_volume == NULL )
    {
        kvsMessageError( "volume is not setted for the thread" );
        return false;
    }

    if ( m_locator == NULL )
    {
        kvsMessageError( "locator is not setted for the thread" );
        return false;
    }

    return true;
}

void HyperStreamlineThread::run()
{

#ifdef DEBUG
    if ( !check() )
        return;
#endif

    std::vector<float> coords;
    std::vector<unsigned char> colors;

    //calculate the line from the seed point
    const kvs::Vector3f seed_vector = this->calculate_vector( m_seed_point );
    
    switch ( m_direction )
    {
    case HyperStreamline::ForwardDirection:
        {       
            this->calculate_one_side(
                    coords,
                    colors,
                    m_seed_point,
                    seed_vector );
            break;
        }
    case HyperStreamline::BackwardDirection:
        {
            this->calculate_one_side(
            coords,
            colors,
            m_seed_point,
            -seed_vector );
            break;
        }
    case HyperStreamline::BothDirections:
        {
             // Forward direction.
            std::vector<float> tmp_coords1;
            std::vector<unsigned char> tmp_colors1;
            this->calculate_one_side(
                     tmp_coords1,
                     tmp_colors1,
                     m_seed_point,
                     seed_vector ); 

            // backward direction.
            std::vector<kvs::Real32> tmp_coords2;
            std::vector<kvs::UInt8> tmp_colors2;
            this->calculate_one_side(
                    tmp_coords2,
                    tmp_colors2,
                    m_seed_point,
                    -seed_vector ); 

#ifdef MISES_COLOR
            const size_t nvertices1 = tmp_coords1.size() / 3;
            for( size_t i = 0; i < nvertices1; i++ )
            {
                const size_t id  = nvertices1 - i - 1;
                const size_t id3 = 3 * id;

                coords.push_back( tmp_coords1[id3  ] );
                coords.push_back( tmp_coords1[id3+1] );
                coords.push_back( tmp_coords1[id3+2] );
                colors.push_back( tmp_colors1[id3]   );
                colors.push_back( tmp_colors1[id3+1] );
                colors.push_back( tmp_colors1[id3+2] );
            }

            const size_t nvertices2 = tmp_coords2.size() / 3;
            for( size_t i = 1; i < nvertices2; i++ )
            {
                const size_t id3 = 3 * i;

                coords.push_back( tmp_coords2[id3  ] );
                coords.push_back( tmp_coords2[id3+1] );
                coords.push_back( tmp_coords2[id3+2] );
                colors.push_back( tmp_colors2[id3]   );
                colors.push_back( tmp_colors2[id3+1] );
                colors.push_back( tmp_colors2[id3+2] );
            }
#else
            const size_t nvertices1 = tmp_coords1.size() / 3;
            for( size_t i = 0; i < nvertices1; i++ )
            {
                const size_t id  = nvertices1 - i - 1;
                const size_t id3 = 3 * id;

                coords.push_back( tmp_coords1[id3  ] );
                coords.push_back( tmp_coords1[id3+1] );
                coords.push_back( tmp_coords1[id3+2] );

                if ( i < nvertices1 / 2 )
                {
                    float temp = m_eigenvalues[i];
                    m_eigenvalues[i] = m_eigenvalues[id];
                    m_eigenvalues[id] = temp;
                }
            }

            const size_t nvertices2 = tmp_coords2.size() / 3;
            for( size_t i = 1; i < nvertices2; i++ )
            {
                const size_t id3 = 3 * i;

                coords.push_back( tmp_coords2[id3  ] );
                coords.push_back( tmp_coords2[id3+1] );
                coords.push_back( tmp_coords2[id3+2] );
            }
#endif
            break;
        }
    default:
        {
            break;
        }
    }

    m_line.setCoords( kvs::ValueArray<float>(coords) );
    //m_line.setColors( kvs::ValueArray<unsigned char>(colors) );
    m_line.setColorType( kvs::LineObject::VertexColor );
    m_line.setLineType( kvs::LineObject::Uniline );

}

const bool HyperStreamlineThread::calculate_one_side(
    std::vector<kvs::Real32>& coords,
    std::vector<kvs::UInt8>& colors,
    const kvs::Vector3f& seed_point,
    const kvs::Vector3f& seed_vector )
{
    // Register the seed point.
    kvs::Vector3f current_vertex = seed_point;
    kvs::Vector3f next_vertex = seed_point;

    coords.push_back( seed_point.x() );
    coords.push_back( seed_point.y() );
    coords.push_back( seed_point.z() );

    // Register the vector on the seed point.
    kvs::Vector3f current_vector = seed_vector;
    kvs::Vector3f previous_vector = seed_vector;

    // Set the color of seed point.
    kvs::RGBColor col = this->calculate_color( );

    colors.push_back( col.r() );
    colors.push_back( col.g() );
    colors.push_back( col.b() );

    size_t integral_times = 0;
    for ( ; ; )
    {
        // Calculate the next vertex.
        if ( !this->calculate_next_vertex(
                 current_vertex,
                 current_vector,
                 &next_vertex ) )
        {
            return( true );
        }

        // Check the termination.
        if ( this->check_for_termination(
                 current_vertex,
                 current_vector,
                 integral_times,
                 next_vertex ) )
        {
            return( true );
        }

        // Update the vertex and vector.
        current_vertex = next_vertex;
        previous_vector = current_vector;

        coords.push_back( current_vertex.x() );
        coords.push_back( current_vertex.y() );
        coords.push_back( current_vertex.z() );

        // Interpolate vector from vertex of cell.
        float eigenvalue = 0.0f;
        current_vector = this->interpolate_vector( current_vertex, previous_vector, eigenvalue );
        m_eigenvalues.push_back( eigenvalue );

        // Set color of vertex.
        kvs::RGBColor col = this->calculate_color( );

        colors.push_back( col.r() );
        colors.push_back( col.g() );
        colors.push_back( col.b() );

        integral_times++;
    }
}

const bool HyperStreamlineThread::calculate_next_vertex(
    const kvs::Vector3f& current_vertex,
    const kvs::Vector3f& current_direction,
    kvs::Vector3f* next_vertex)
{
    switch( m_method )
    {
    case HyperStreamline::Euler:
        {
            return( this->integrate_by_euler(
                    current_vertex,
                    current_direction,
                    &(*next_vertex) ) );
            break;
        }
    case HyperStreamline::RungeKutta2nd:
        {
            return( this->integrate_by_runge_kutta_2nd(
                    current_vertex,
                    current_direction,
                    &(*next_vertex) ) );
            break;
        }
    case HyperStreamline::RungeKutta4th:
        {
            return( this->integrate_by_runge_kutta_4th(
                    current_vertex,
                    current_direction,
                    &(*next_vertex) ) );
            break;
        }
    default: break;
    }

    kvsMessageError( "Specified integral method is not defined." );
    return( false );
}

const bool HyperStreamlineThread::integrate_by_euler(
    const kvs::Vector3f& current_vertex,
    const kvs::Vector3f& current_direction,
    kvs::Vector3f* next_vertex )
{
    if ( m_enable_boundary_condition )
    {
        if ( !this->check_for_inside_volume( current_vertex ) ) return( false );
    }

    const kvs::Vector3f k1 = m_integration_interval * current_direction;
    *next_vertex = current_vertex + k1;

    return( true );
}

const bool HyperStreamlineThread::integrate_by_runge_kutta_2nd(
    const kvs::Vector3f& current_vertex,
    const kvs::Vector3f& current_direction,
    kvs::Vector3f* next_vertex )
{
    if ( m_enable_boundary_condition )
    {
        if ( !this->check_for_inside_volume( current_vertex ) ) return( false );
    }

    const kvs::Vector3f k1 = m_integration_interval * current_direction;
    // Interpolate vector from vertex of cell.
    const kvs::Vector3f vertex = current_vertex + 0.5f * k1;

    if ( m_enable_boundary_condition )
    {
        if ( !this->check_for_inside_volume( vertex ) ) return( false );
    }

    float egv = 0.0f;
    const kvs::Vector3f direction = this->interpolate_vector( vertex, current_direction, egv );
    const kvs::Vector3f k2 = m_integration_interval * direction;
    *next_vertex = vertex + k2;

    return( true );
}

const bool HyperStreamlineThread::integrate_by_runge_kutta_4th(
    const kvs::Vector3f& current_vertex,
    const kvs::Vector3f& current_direction,
    kvs::Vector3f* next_vertex )
{
    if ( m_enable_boundary_condition )
    {
        if ( !this->check_for_inside_volume( current_vertex ) ) return( false );
    }

    // Calculate integration interval.
    const float interval = m_integration_interval / static_cast<float>(current_direction.length());

    const kvs::Vector3f k1 = interval * current_direction;

    // Interpolate vector from vertex of cell.
    const kvs::Vector3f vertex2 = current_vertex + 0.5f * k1;

    if ( m_enable_boundary_condition )
    {
        if ( !this->check_for_inside_volume( vertex2 ) ) return( false );
    }

    float egv = 0.0f;
    const kvs::Vector3f direction2 = this->interpolate_vector( vertex2, current_direction, egv );
    const kvs::Vector3f k2 = interval * direction2;

    // Interpolate vector from vertex of cell.
    const kvs::Vector3f vertex3 = current_vertex + 0.5f * k2;

    if ( m_enable_boundary_condition )
    {
        if ( !this->check_for_inside_volume( vertex3 ) ) return( false );
    }

    const kvs::Vector3f direction3 = this->interpolate_vector( vertex3, current_direction, egv );
    const kvs::Vector3f k3 = interval * direction3;

    // Interpolate vector from vertex of cell.
    const kvs::Vector3f vertex4 = current_vertex + k3;

    if ( m_enable_boundary_condition )
    {
        if ( !this->check_for_inside_volume( vertex4 ) ) return( false );
    }

    const kvs::Vector3f direction4 = this->interpolate_vector( vertex4, current_direction, egv );
    const kvs::Vector3f k4 = interval * direction4;

    *next_vertex = current_vertex + ( k1 + 2.0f * ( k2 + k3 ) + k4 ) / 6.0f;

    return( true );
}


const bool HyperStreamlineThread::check_for_inside_volume( const kvs::Vector3f& point )
{
    switch( m_volume->volumeType() )
    {
        case kvs::VolumeObjectBase::Structured:
        {
            const kvs::StructuredVolumeObject* structured_volume =
            reinterpret_cast<const kvs::StructuredVolumeObject*>( m_volume );

            const float dimx = static_cast<float>( structured_volume->resolution().x() - 1 );
            const float dimy = static_cast<float>( structured_volume->resolution().y() - 1 );
            const float dimz = static_cast<float>( structured_volume->resolution().z() - 1 );

            if ( point.x() < 0.0f || dimx < point.x() ) return( false );
            if ( point.y() < 0.0f || dimy < point.y() ) return( false );
            if ( point.z() < 0.0f || dimz < point.z() ) return( false );

            return( true );
        }
        case kvs::VolumeObjectBase::Unstructured:
        {
            if ( point.x() < m_volume->minExternalCoord().x() || point.x() > m_volume->maxExternalCoord().x() ) return false;
            if ( point.y() < m_volume->minExternalCoord().y() || point.y() > m_volume->maxExternalCoord().y() ) return false;
            if ( point.z() < m_volume->minExternalCoord().z() || point.z() > m_volume->maxExternalCoord().z() ) return false;

            float pos[3] = { point.x(), point.y(), point.z() };
            if ( m_locator->findCell( pos ) != -1 )
                return true;
            else
                return false;
        }
    }

    return false;

}

const bool HyperStreamlineThread::check_for_vector_length( const kvs::Vector3f& direction )
{
    return( direction.length() < m_vector_length_threshold );
}

const bool HyperStreamlineThread::check_for_integration_times( const size_t times )
{
    return( times <= m_integration_times_threshold );
}


const bool HyperStreamlineThread::check_for_acceptance( const std::vector<kvs::Real32>& vertices )
{
    kvs::IgnoreUnusedVariable( vertices );
    return( true );
}

const bool HyperStreamlineThread::check_for_termination(
    const kvs::Vector3f& current_vertex,
    const kvs::Vector3f& direction,
    const size_t         integration_times,
    const kvs::Vector3f& next_vertex )
{
    kvs::IgnoreUnusedVariable( current_vertex );

    if ( m_enable_boundary_condition )
    {
        if ( !check_for_inside_volume( next_vertex ) ) 
            return( true );
    }

    if ( m_enable_vector_length_condition )
    {
        if ( !check_for_vector_length( direction ) )
            return true;
    }

    if ( m_enable_integration_times_condition )
    {
        if ( !check_for_integration_times( integration_times ) )
            return true;
    }

    //if ( m_degenerate )
    //{
    //    return true;
    //}

    return( false );
}

const kvs::Vector3f HyperStreamlineThread::calculate_vector( const kvs::Vector3f& point )
{
    const kvs::Vector3f origin( 0.0f, 0.0f, 0.0f );
    float egv = 0.0f;
    const kvs::Vector3f seed_vector = this->interpolate_vector( point, origin, egv );
    m_eigenvalues.push_back( egv );
    return( seed_vector );
}


const kvs::RGBColor HyperStreamlineThread::calculate_color()
{
#ifdef MISES_COLOR
    return( m_transfer_function.colorMap().at(m_mises) );
#else
    return( kvs::RGBColor( 255, 0, 0 ) );
#endif
}


const kvs::Vector3f HyperStreamlineThread::interpolate_vector(
    const kvs::Vector3f& vertex,
    const kvs::Vector3f& previous_vector,
    float& eigenvalue )
{
    //tensor field interpolation for unstructured grid
    if ( m_volume->volumeType() == kvs::VolumeObjectBase::Unstructured )   
    {
        const kvs::UnstructuredVolumeObject* volume = static_cast<const kvs::UnstructuredVolumeObject*>( m_volume );

        // locate the cell that contains the vertex
        float gpos[3] = { vertex[0], vertex[1], vertex[2] };

        int cellid = m_locator->findCell( gpos ); 

#ifdef DEBUG
        if ( cellid < 0 || cellid > volume->connections().size() )
        {
            std::cerr << "Wrong, cellid is not right!!" << std::endl;
        }
#endif

        // interpolate the tensor matrix 
        int length = m_volume->cellType();
        const float *intfunc = m_locator->cell()->interpolationFunctions( m_locator->cell()->localPoint() );
        KVS_ASSERT( cellid != -1 );
        kvs::Matrix33f stress_interpolated; // initilaze to zero

        m_mises = 0;
        kvs::Matrix33f stress_each_vertex;

        if ( m_volume->veclen() == 7 )
        {
             for ( int i = 0 ; i < length; i ++ )
            {
                int vi = volume->connections()[ cellid * length + i ];  //vertex index
                const float* s = static_cast<const float*>( volume->values().pointer() );         

                stress_each_vertex = kvs::Matrix33f( 
                    s[ vi * 7 ],     s[ vi * 7 + 3 ], s[ vi * 7 + 4 ],
                    s[ vi * 7 + 3 ], s[ vi * 7 + 1 ], s[ vi * 7 + 5 ],
                    s[ vi * 7 + 4 ], s[ vi * 7 + 5 ], s[ vi * 7 + 2 ] ); 

                stress_each_vertex *= intfunc[i];
                stress_interpolated += stress_each_vertex;  
                m_mises += intfunc[i] * s[ vi * 7 + 6 ];
            }
        }
        else if ( m_volume->veclen() == 6 )
        {
            for ( int i = 0 ; i < length; i ++ )
            {
                int vi = volume->connections()[ cellid * length + i ];  //vertex index
                const float* s = static_cast<const float*>( volume->values().pointer() );         

                stress_each_vertex = kvs::Matrix33f( 
                    s[ vi * 6 ],     s[ vi * 6 + 3 ], s[ vi * 6 + 4 ],
                    s[ vi * 6 + 3 ], s[ vi * 6 + 1 ], s[ vi * 6 + 5 ],
                    s[ vi * 6 + 4 ], s[ vi * 6 + 5 ], s[ vi * 6 + 2 ] ); 

                stress_each_vertex *= intfunc[i];
                stress_interpolated += stress_each_vertex;  
            }

#ifdef EIGEN_VALUE_COLOR
            m_mises = 0.0;
#else
            const kvs::Matrix33f& m = stress_each_vertex;
            m_mises = std::sqrt( 
                0.5 * ( 
                std::pow( ( m[0][0]-m[1][1] ), 2 ) + std::pow( ( m[1][1]-m[2][2] ), 2 )
                + std::pow( ( m[2][2]-m[0][0] ), 2 ) + 3 * ( 2*m[0][1]*m[0][1] + 2*m[0][2]*m[0][2] + 
                2*m[1][2]*m[1][2] ) 
                ) 
                );
#endif
                
        }



        if ( kvs::Math::Abs(stress_interpolated.determinant()) < 0.0001f ) // degenerate case
            m_degenerate = true;        
        

        // singular value Decomposition to get the biggest vector
        kvs::EigenDecomposer<float> ed( stress_interpolated, kvs::EigenDecomposer<float>::Symmetric );

        /*
        const kvs::Vector<float>& egvalue = ed.eigenValues();
        float max_index = -1;
        float max = 0;
        for ( int i = 0; i < 3; i ++ )
        {
            if ( kvs::Math::Abs<float>( egvalue[i] ) > max )
            {
                max = kvs::Math::Abs<float>( egvalue[i] );
                max_index = i;
            }
        }
        kvs::Vector3f current_vector( ed.eigenVector( max_index )[0], ed.eigenVector( max_index )[1], ed.eigenVector( max_index )[2] );
        */
        eigenvalue = ed.eigenValue(m_nth_egvector);
        kvs::Vector3f current_vector( ed.eigenVector(m_nth_egvector)[0], ed.eigenVector(m_nth_egvector)[1], ed.eigenVector(m_nth_egvector)[2] );
        
        if ( previous_vector != kvs::Vector3f( 0, 0, 0 ) )
        {
            float cos = current_vector.dot( previous_vector ) / ( current_vector.length() * previous_vector.length() );
            if ( cos > 0.1 )
                return current_vector;
            else if ( cos < -0.1 )
                return -current_vector;
            else
            {
                m_degenerate = true;
                return current_vector;
            }
        }
        else // the seed vector case
        {
            kvs::Vector3f temp_vector( 1, 1, 1 );
            float cos = current_vector.dot( temp_vector ) / ( current_vector.length() * temp_vector.length() );
            if ( cos > 0.1 )
                return current_vector;
            else if ( cos < -0.1 )
                return -current_vector;
            else
            {
                m_degenerate = true;
                return current_vector;
            }
        }
    }
    else
    {
        std::cerr << "not supported yet " << std::endl;
    }
}

void HyperStreamlineThread::setVolume( const kvs::VolumeObjectBase* volume )
{
    m_volume = volume;
}

void HyperStreamlineThread::setTransferFunction( const kvs::TransferFunction& transfer_function )
{
    m_transfer_function = transfer_function;
}

void HyperStreamlineThread::setSeedPoint( const kvs::Vector3f seed )
{
    m_seed_point = seed;
}

void HyperStreamlineThread::setLocator( const kvs::CellTree* ct, const kvs::UnstructuredVolumeObject* volume )
{
    m_locator = new kvs::CellLocatorBIH();
    m_locator->setCellTree( const_cast<kvs::CellTree*>(ct) );
    m_locator->setDataSet( volume );
    m_locator->setMode( kvs::CellLocator::CACHEFULL );
    m_locator->initializeCell();
}

void HyperStreamlineThread::setDirection( const kvs::HyperStreamline::IntegrationDirection Direction)
{
    m_direction = Direction;  
}

void HyperStreamlineThread::setMethod( const kvs::HyperStreamline::IntegrationMethod Method )
{
    m_method = Method;
}

const kvs::LineObject& HyperStreamlineThread::line() const
{
    return m_line;
}

const std::vector<float>& HyperStreamlineThread::eigenValue() const
{
    return m_eigenvalues;
}

void HyperStreamlineThread::setGoWithNthEigenVector( const int nth )
{
    if ( nth < 0 || nth > 2 )
    {
        std::cout << "From HyperStreamline::setGoWithNthEigenVector: There are only 3 eigen vectors!! set to default(0) " << std::endl;
        m_nth_egvector = 0;
        return;
    }
    m_nth_egvector = nth;
}

void HyperStreamlineThread::setDisableCache()
{
    m_locator->setMode( kvs::CellLocator::CACHEOFF );
}

} // end of namespace kvs
