/*
 *  MainScreen.cpp
 *  
 *
 *  Created by njun on 11/10/11.
 *  Copyright 2011 Jun Nishimura. All rights reserved.
 *
 */

#include "MainScreen.h"

#include <kvs/PointRenderer>
#include <kvs/Bounds>

//#define UNIFORM_DISTRIBUTION

bool flag_background=true;

unsigned int opacity_polygon = 10;
unsigned int opacity_line    = 60;
unsigned int opacity_point   = 40;

namespace kvs
{

class Argument : public kvs::CommandLine
{

public:

    Argument( int argc, char** argv ) :
        CommandLine( argc, argv )
    {
        add_help_option();
        add_option( "r", "[size_t] Repeat level. ( default : 1 )", 1, false );

        add_option( "volume1", "[string] kvs::UnstructuredVolumeObject file path. ( optional )", 1, false );
        add_option( "volume2", "[string] kvs::UnstructuredVolumeObject file path. ( optional )", 1, false );
        add_option( "polygon", "[string] kvs::PolygonObject file path. ( optional )", 1, false );
        add_option( "tfunc", "[string] kvs::TransferFunction file path. ( optional )", 1, false );

        add_option( "DisableShading", "Disable shading. ( default : eable shading )", 0, false );

        if( !this->parse() ) exit( EXIT_FAILURE );
    }

};


MainScreen::MainScreen( kvs::glut::Application* app ) :
    kvs::glut::Screen( app )
{
    this->initialize( app );
}

MainScreen::~MainScreen( void )
{
    this->clear();
}

void MainScreen::redraw()
{
    if ( flag_background )
        this->background()->setColor( kvs::RGBAColor( 255,255,255,1 ) );
    else
        this->background()->setColor( kvs::RGBAColor( 0,0,0,1 ) );

    BaseClass::redraw();
}

void MainScreen::initialize( kvs::glut::Application* app )
{
    Argument arg( app->argc(), app->argv() );
    BaseClass::setTitle( "kvs::MainScreen" );
    BaseClass::setGeometry( 0, 0, 512, 512 );
    
    //volume1
    std::string default_volume1( "../../data/engine/v6engine_stress_and_mises.kvsml" );
    if ( !arg.hasOption( "volume1" ) )
    {
        kvsMessageWarning( m_volume1, "no volume inputed, default volume used");
        m_volume1 = new kvs::UnstructuredVolumeImporter( default_volume1 );
    }
    else
        m_volume1 = new kvs::UnstructuredVolumeImporter( arg.optionValue<std::string>( "volume1" ) ); 
    if ( !m_volume1 )
    {
        kvsMessageError( "Cannot create an unstructured volume object." );
        exit( EXIT_FAILURE );
    }

    //volume2
    std::string default_volume2( "../../data/engine/v6engine.kvsml" );
    if ( !arg.hasOption( "volume2" ) )
    {
        kvsMessageWarning( m_volume2, "no volume inputed, default volume used");
        m_volume2 = new kvs::UnstructuredVolumeImporter( default_volume2 );
    }
    else
        m_volume2 = new kvs::UnstructuredVolumeImporter( arg.optionValue<std::string>( "volume2" ) ); 
    if ( !m_volume2 )
    {
        kvsMessageError( "Cannot create an unstructured volume object." );
        exit( EXIT_FAILURE );
    }

    //cell locator
    m_locator = new kvs::CellLocatorBIH();
    m_locator->setDataSet( m_volume1 );
    m_locator->setMode( kvs::CellLocator::CACHEHALF );
    m_locator->setParallel();
    m_locator->initializeCell();
    m_locator->build();

    //polygon
    if ( !arg.hasOption( "polygon" ) )
    {
        kvsMessageWarning( m_polygon, "no polygon inputed, generate");
        kvs::PolygonObject* import_polygon = new kvs::ExternalFaces( m_volume2 );
        m_polygon = new kvs::PolygonToPolygon( import_polygon );
    }
    else
    {
        kvs::PolygonObject* import_polygon = new kvs::PolygonImporter( arg.optionValue<std::string>( "polygon" ) );
        m_polygon = new kvs::PolygonToPolygon( import_polygon );
    }
    if ( !m_polygon )
    {
        kvsMessageError( "Cannot create surface." );
        exit( EXIT_FAILURE );
    }

#ifdef WRITE
    kvs::KVSMLObjectPolygon* kvsml_polygon = new kvs::PolygonExporter<kvs::KVSMLObjectPolygon>( m_polygon );
    kvsml_polygon->setWritingDataType( kvs::KVSMLObjectPolygon::ExternalBinary );
    kvsml_polygon->write( "pump_external.kvsml" );
#endif

    //seed points
    float xmax = m_volume1->maxExternalCoord().x();
    float xmin = m_volume1->minExternalCoord().x();
    float ymax = m_volume1->maxExternalCoord().y();
    float ymin = m_volume1->minExternalCoord().y();
    float zmax = m_volume1->maxExternalCoord().z();
    float zmin = m_volume1->minExternalCoord().z();

#ifdef UNIFORM_DISTRIBUTION
    unsigned int xres = 20;
    unsigned int yres = 20;
    unsigned int zres = 20;

    std::vector<kvs::Real32>v;
    for( unsigned int i = 0; i < xres; i++ )
    {
        for ( unsigned int j = 0; j < yres; j ++ )
        {
            for ( unsigned int k = 0; k < zres; k ++ )
            {
                v.push_back( (xmax-xmin)/xres*i + xmin );
                v.push_back( (ymax-ymin)/yres/2*j + ymin );
                v.push_back( zmax - (zmax-zmin)/zres/2*k );
            }
        }   
    }
    m_seed_point = new kvs::PointObject( kvs::ValueArray<kvs::Real32>( v ) );
#else
    const kvs::Vector3f center = m_volume1->objectCenter();
    unsigned int xres = 100;
    unsigned int yres = 100;
    unsigned int zres = 100;

    std::vector<kvs::Real32>v;
    for( int i = -1; i < 2; i++ )
    {
        for ( int j = -1; j < 2; j ++ )
        {
            for ( int k = -1; k < 2; k ++ )
            {
                v.push_back( (xmax-xmin)/xres*i + center.x() );
                v.push_back( (ymax-ymin)/yres*j + center.y() );
                v.push_back( (zmax-zmin)/zres*k + center.z() );
            }
        }   
    }
    
    kvs::ValueArray<kvs::Real32>coords(3);
    coords[0] = center.x();
    coords[1] = center.y();
    coords[2] = center.z();
    m_seed_point = new kvs::PointObject( coords );
#endif
 
    m_seed_point->updateMinMaxCoords();

    if ( !m_seed_point )
    {
        kvsMessageError( "Cannot creat a point object." );
        delete m_volume1;
        exit( EXIT_FAILURE );
    }
    m_seed_point->setColor( kvs::RGBColor( 255, 0, 0 ) );
    m_seed_point->setSize( 10.0f );
    m_seed_point->setMinMaxObjectCoords( m_volume1->minObjectCoord(), m_volume1->maxObjectCoord() );
    m_seed_point->setMinMaxExternalCoords( m_volume1->minExternalCoord(), m_volume1->maxExternalCoord() );

    m_tfunc.adjustRange( m_volume1 );
    m_streamline = new kvs::HyperStreamline();
	m_streamline->setGoWithNthEigenVector( 0 );
	m_streamline->setIntegrationDirection( kvs::HyperStreamline::BothDirections );
	m_streamline->setLocator( m_locator->cellTree(), m_volume1 );
	m_streamline->setIntegrationMethod( kvs::HyperStreamline::RungeKutta2nd );
    m_streamline->setTransferFunction( m_tfunc );
	m_streamline->setSeedPoints( m_seed_point );
	m_streamline->exec( m_volume1 );
    m_buffered_streamline = new kvs::HyperStreamline();
    m_buffered_streamline->setTransferFunction( m_tfunc );
    m_buffered_streamline->setColors( m_streamline->colors() );
    m_buffered_streamline->setCoords( m_streamline->coords() );
    m_buffered_streamline->setConnections( kvs::ValueArray<kvs::UInt32>(0) );

    if ( !m_streamline )
    {
        kvsMessageError( "Cannot creat a streamline object." );
        delete m_volume1;
        delete m_seed_point;
        exit( EXIT_FAILURE );
    }
    m_streamline->setName( "Streamline" );

#ifdef WRITE
    kvs::KVSMLObjectLine* kvsml_line = new kvs::LineExporter<kvs::KVSMLObjectLine>( m_streamline );
    kvsml_line->setWritingDataType( kvs::KVSMLObjectLine::ExternalBinary );
    kvsml_line->write( "pump_line.kvsml" );
#endif

    //renderer
    const size_t repeat_level = arg.hasOption( "r" ) ? arg.optionValue<size_t>( "r" ) : 1;
    m_renderer = new kvs::glew::StochasticRenderer( repeat_level );
    m_renderer->enableLODControl();

    kvs::NullObject* null = new::kvs::NullObject( m_volume2 );

    m_volume_renderer = new kvs::glew::StochasticVolumeRenderer( m_volume2 );
    m_volume_renderer->setShader( kvs::Shader::BlinnPhong() );
    m_renderer->registerRenderer( m_volume_renderer );

    m_polygon->setOpacity( opacity_polygon );
    m_polygon_renderer = new kvs::glew::StochasticPolygonRenderer( m_polygon );
    m_polygon_renderer->setShader( kvs::Shader::BlinnPhong() );
    m_renderer->registerRenderer( m_polygon_renderer );

    m_point_renderer = new kvs::glew::StochasticPointRenderer( m_seed_point );
#ifndef UNIFORM_DISTRIBUTION
    m_renderer->registerRenderer( m_point_renderer );
#endif

    m_line_renderer = new kvs::glew::StochasticLineRenderer( m_streamline );
    m_line_renderer->setShader( kvs::Shader::BlinnPhong() );
    m_line_renderer->setOpacity( opacity_line );
    m_renderer->registerRenderer( m_line_renderer );

    m_buffered_line_renderer = new kvs::glew::StochasticLineRenderer( m_buffered_streamline );
    m_buffered_line_renderer->setShader( kvs::Shader::BlinnPhong() );
    m_buffered_line_renderer->setOpacity( opacity_line );
    m_renderer->registerRenderer( m_buffered_line_renderer );

    BaseClass::registerObject( null, m_renderer );

}

void MainScreen::keyPressEvent( kvs::KeyEvent* event )
{
    BaseClass::keyPressEvent( event );
    
    switch ( event->key() )
    {
        case kvs::Key::o: this->controlTarget() = kvs::ScreenBase::TargetObject; break;
        case kvs::Key::l: this->controlTarget() = kvs::ScreenBase::TargetLight; break;
        case kvs::Key::b: flag_background = ~flag_background; break;
        case kvs::Key::Up:
        {
            if ( opacity_polygon < 255 )
                opacity_polygon += 5;
            m_polygon->setOpacity( opacity_polygon );
            m_renderer->clearEnsembleBuffer();
            this->redraw();
            break;
        }
        case kvs::Key::Down:
        {
            if ( opacity_polygon > 0 )
                opacity_polygon -= 5;
            m_polygon->setOpacity( opacity_polygon );
            m_renderer->clearEnsembleBuffer();
            this->redraw();
            break;
        }
        case kvs::Key::Left:
        {
            if ( opacity_line < 255 )
                opacity_line += 5;
            m_line_renderer->setOpacity( opacity_line );
            m_renderer->clearEnsembleBuffer();
            this->redraw();
            break;
        }
        case kvs::Key::Right:
        {
            if ( opacity_line > 0 )
                opacity_line -= 5;
            m_line_renderer->setOpacity( opacity_line );
            m_renderer->clearEnsembleBuffer();
            this->redraw();
            break;
        }
        case kvs::Key::s:
        {
            m_buffer.push_back( m_streamline->connections()[1] );

            const kvs::ValueArray<kvs::UInt8>& temp_colors = m_streamline->colors();
            const kvs::ValueArray<kvs::UInt32>& temp_connections = m_streamline->connections();
            const kvs::ValueArray<kvs::Real32>& temp_coords = m_streamline->coords();

            const kvs::ValueArray<kvs::UInt8>& buf_colors = m_buffered_streamline->colors();
            const kvs::ValueArray<kvs::UInt32>& buf_connections = m_buffered_streamline->connections();
            const kvs::ValueArray<kvs::Real32>& buf_coords = m_buffered_streamline->coords();

            const unsigned int size_color = temp_colors.size() + buf_colors.size();
            const unsigned int size_conne = temp_connections.size() + buf_connections.size();
            const unsigned int size_coord = temp_coords.size() + buf_coords.size();

            kvs::ValueArray<kvs::UInt8> colors( size_color );
            kvs::ValueArray<kvs::UInt32> connections( size_conne );
            kvs::ValueArray<kvs::Real32> coords( size_coord );

            for ( unsigned int i = 0; i < temp_colors.size(); i ++ )
                colors[i] = buf_colors[i];
            for ( unsigned int i = temp_colors.size(); i < size_color; i ++ )
                colors[i] = temp_colors[i];

            for ( unsigned int i = 0; i < temp_connections.size(); i ++ )
                connections[i] = buf_connections[i];
            for ( unsigned int i = temp_connections.size(); i < size_conne; i ++ )
                connections[i] = temp_connections[i];

            for ( unsigned int i = 0; i < temp_coords.size(); i ++ )
                coords[i] = buf_coords[i];
            for ( unsigned int i = temp_coords.size(); i < size_coord; i ++ )
                coords[i] = temp_coords[i];

            kvs::HyperStreamline* new_line = new kvs::HyperStreamline();
            new_line->setColors( colors );
            new_line->setConnections( connections );
            new_line->setCoords( coords );

            m_renderer->changeObject( new_line, m_buffered_line_renderer, true );
            m_buffered_streamline = new_line;
            break;
        }
        case kvs::Key::Backslash:
        {
            if ( m_buffered_streamline->connections().size() > 1 )
            {
                const unsigned int interval = m_buffer.back();
                m_buffer.pop_back();

                kvs::ValueArray<kvs::UInt8> buf_colors;
                buf_colors.deepCopy(m_buffered_streamline->colors());
                kvs::ValueArray<kvs::UInt32> buf_connections;
                buf_connections.deepCopy(m_buffered_streamline->connections());
                kvs::ValueArray<kvs::Real32> buf_coords;
                buf_coords.deepCopy(m_buffered_streamline->coords());
                
                const unsigned int size_color = buf_colors.size()-interval*3;
                const unsigned int size_conne = buf_connections.size()-2;
                const unsigned int size_coord = buf_coords.size()-interval*3;

                kvs::ValueArray<kvs::UInt8> colors( size_color );
                kvs::ValueArray<kvs::UInt32> connections( size_conne );
                kvs::ValueArray<kvs::Real32> coords( size_coord );

                for ( unsigned int i = 0; i < size_color; i ++ )
                    colors[i] = buf_colors[i];

                for ( unsigned int i = 0; i < size_conne; i ++ )
                    connections[i] = buf_connections[i];

                for ( unsigned int i = 0; i < size_coord; i ++ )
                    coords[i] = buf_coords[i];

                kvs::HyperStreamline* new_line = new kvs::HyperStreamline();
                new_line->setColors( colors );
                new_line->setConnections( connections );
                new_line->setCoords( coords );

                m_renderer->changeObject( new_line, m_buffered_line_renderer, true );
                m_buffered_streamline = new_line;
            }
            break;
        }
        default: break;
    }
}

void MainScreen::clear( void )
{
    if ( m_volume1 ) delete m_volume1;
}

kvs::PointObject* MainScreen::seedPoint( void )
{
    return( m_seed_point );
}

kvs::UnstructuredVolumeObject* MainScreen::volume1()
{
    return( m_volume1 );
}

kvs::UnstructuredVolumeObject* MainScreen::volume2()
{
    return( m_volume2 );
}

kvs::glew::StochasticVolumeRenderer* MainScreen::volumeRenderer()
{
    return m_volume_renderer;
}

kvs::glew::StochasticRenderer* MainScreen::renderer()
{
    return m_renderer;
}

void MainScreen::updateStreamLine( void )
{
    kvs::HyperStreamline* streamline = new kvs::HyperStreamline();
	streamline->setGoWithNthEigenVector( 0 );
	streamline->setIntegrationDirection( kvs::HyperStreamline::BothDirections );
	streamline->setLocator( m_locator->cellTree(), m_volume1 );
	streamline->setIntegrationMethod( kvs::HyperStreamline::RungeKutta2nd );
    streamline->setTransferFunction( m_tfunc );
	streamline->setSeedPoints( m_seed_point );
	streamline->exec( m_volume1 );

    if ( !m_streamline )
    {
        kvsMessageError( "Cannot create a streamline object." );
        delete m_volume1;
        delete m_seed_point;
        exit( EXIT_FAILURE );
    }

    m_renderer->changeObject( streamline, m_line_renderer, true );
    m_renderer->changeObject( m_seed_point, m_point_renderer, false );

    m_streamline = streamline;

}

}
