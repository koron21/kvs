#include <kvs/glut/Application>
#include <kvs/glut/Screen>

#include <kvs/CommandLine>
#include <kvs/ExternalFaces>
#include <kvs/LineImporter>
#include <kvs/PolygonImporter>
#include <kvs/RGBFormulae>
#include <kvs/TimerEventListener>
#include <kvs/UnstructuredVolumeImporter>

#include "NullObject.h"
#include "PolygonToPolygon.h"
#include "StochasticLineRenderer.h"
#include "StochasticPointRenderer.h"
#include "StochasticPolygonRenderer.h"
#include "StochasticRenderer.h"
#include "StochasticVolumeRenderer.h"

kvs::LineObject* line = NULL;
kvs::PointObject* point = NULL;
kvs::glew::StochasticRenderer* renderer = NULL;
kvs::glew::StochasticPointRenderer* point_renderer= NULL;
float speed = 0.5;
float time_step = 0;
unsigned int max_length = 0;

bool stop_flag = false;
bool restart_flag = false;

class Argument : public kvs::CommandLine
{

public:

    Argument( int argc, char** argv ) :
        CommandLine( argc, argv )
    {
        add_option( "line", "[string] kvs::PolygonObject file path. ( optional )", 1, true );
        add_option( "volume", "[string] kvs::UnstructuredVolumeObject file path. ( optional )", 1, false );
        add_option( "polygon", "[string] kvs::PolygonObject file path. ( optional )", 1, false );
        add_option( "speed", "[float] specify tracing speed multiplier. ( optional )", 1, 1.0f );
        add_help_option();
        if( !this->parse() ) exit( EXIT_FAILURE );
    }

};

class KeyPressEvent : public kvs::KeyPressEventListener
{
    void update( kvs::KeyEvent* event )
    {
       if ( event->key() == kvs::Key::s ) 
           stop_flag = !stop_flag;

       if ( event->key() == kvs::Key::d ) 
       {
           time_step = 0;
           renderer->clearEnsembleBuffer();
       }

       if ( event->key() == kvs::Key::c ) 
       {
           renderer->clearEnsembleBuffer();
       }
    }
};

class TimerEvent : public kvs::TimerEventListener
{
public:
    void update( kvs::TimeEvent* event )
    {
        if ( !stop_flag )
        {
            // do some stuff here to generate the new points
            if ( (unsigned int)(time_step) <= max_length )
            {
                const kvs::ValueArray<float>& old_coords = point->coords();
                const kvs::ValueArray<unsigned char>& old_colors = point->colors();
                const unsigned int nlines = line->nconnections();
                const unsigned int nvertices = line->nvertices();
                switch ( line->lineType() )
                {
                    case kvs::LineObject::Polyline:
                    {
                        std::vector<float>coords( nvertices*3 );
                        std::vector<unsigned char>colors( nvertices*3 );

                        for ( unsigned int i = 0; i < nlines; i ++ )
                        {
                            const kvs::ValueArray<float>& line_coords = line->coords();
                            const kvs::ValueArray<unsigned char>& line_colors = line->colors();
                            const unsigned int start = line->connection( i ).x();
                            const unsigned int end = line->connection( i ).y();
                            const unsigned int current = start + (unsigned int)(time_step) % ( end-start+1 );

                            if ( (unsigned int)(time_step) > end - start )
                            {
                                for ( unsigned int j = 0; j <= end-start; j ++ )
                                {
                                    coords.push_back( line_coords[ 3*(start+j) ] );
                                    coords.push_back( line_coords[ 3*(start+j)+1 ] );
                                    coords.push_back( line_coords[ 3*(start+j)+2 ] );

                                    colors.push_back( line_colors[ 3*(start+j) ] );
                                    colors.push_back( line_colors[ 3*(start+j)+1 ] );
                                    colors.push_back( line_colors[ 3*(start+j)+2 ] );
                                }
                            }
                            else
                            {
                                for ( unsigned int j = 0; j <= current-start; j ++ )
                                {
                                    coords.push_back( line_coords[ 3*(start+j) ] );
                                    coords.push_back( line_coords[ 3*(start+j)+1 ] );
                                    coords.push_back( line_coords[ 3*(start+j)+2 ] );

                                    colors.push_back( line_colors[ 3*(start+j) ] );
                                    colors.push_back( line_colors[ 3*(start+j)+1 ] );
                                    colors.push_back( line_colors[ 3*(start+j)+2 ] );
                                }
                            }

                        }

                        const_cast<kvs::ValueArray<float>& >(old_coords).deallocate();
                        const_cast<kvs::ValueArray<unsigned char>& >(old_colors).deallocate();

                        point->setCoords( kvs::ValueArray<float>(coords) );
                        point->setColors( kvs::ValueArray<unsigned char>(colors) );
                        time_step += speed;

                        break;
                    }
                    case kvs::LineObject::Uniline:
                    {
                        break;
                    }
                    default:
                    {
                        std::cout << "Line type not support by filter yet" << std::endl;
                        return;
                    }
                }

                point_renderer->attachObject( point );
                point_renderer->clearEnsembleBuffer();
                //if ( (unsigned int)(time_step/speed ) % 30 == 0 )
                //    renderer->clearEnsembleBuffer();
                point_renderer->enableUpdateFlag();
            }
        }
        screen()->redraw();
    }
};

int main( int argc, char** argv )
{
    // arguments
    kvs::glut::Application app( argc, argv );
    Argument arg( app.argc(), app.argv() );

    // screen
    kvs::glut::Screen screen( &app );
    int interval = 30;
    kvs::glut::Timer timer( interval );
    TimerEvent timer_event;
    KeyPressEvent key_event;
    screen.setTitle( "SPT Particle Tracer" );
    screen.addKeyPressEvent( &key_event );
    screen.addTimerEvent( &timer_event, &timer );
    screen.show();
    screen.showFullScreen();

    // volume
    kvs::UnstructuredVolumeObject* volume = NULL;
    std::string default_volume( "../../data/engine/v6engine_stress_and_mises.kvsml" );

    if ( !arg.hasOption( "volume" ) )
    {
        kvsMessageWarning( volume, "no volume inputed, default volume used");
        volume = new kvs::UnstructuredVolumeImporter( default_volume );
    }
    else
        volume = new kvs::UnstructuredVolumeImporter( arg.optionValue<std::string>( "volume" ) ); 

    // polygon
    kvs::PolygonObject* polygon = NULL;
    if ( !arg.hasOption( "polygon" ) )
    {
        kvsMessageWarning( polygon, "no polygon inputed, generate");
        kvs::PolygonObject* import_polygon = new kvs::ExternalFaces( volume );
        polygon = new kvs::PolygonToPolygon( import_polygon );
    }
    else
    {
        kvs::PolygonObject* import_polygon = new kvs::PolygonImporter( arg.optionValue<std::string>( "polygon" ) );
        polygon = new kvs::PolygonToPolygon( import_polygon );
    }

    // line
    line = new kvs::LineImporter( arg.optionValue<std::string>( "line" ) );
    for ( unsigned int i = 0; i < line->nconnections(); i ++ )
    {
        const unsigned int start = line->connection(i).x();
        const unsigned int end = line->connection(i).y();

        if ( end - start > max_length )
            max_length = end - start;
    }

    // do some stuff here to initialize points
    if ( arg.hasOption( "speed" ) )
        speed = arg.optionValue<float>( "speed" );

    point = new kvs::PointObject();

    const unsigned int nlines = line->nconnections();
    switch ( line->lineType() )
    {
        case kvs::LineObject::Polyline:
        {
            kvs::ValueArray<float>coords( nlines*3 );
            kvs::ValueArray<unsigned char>colors( nlines*3 );

            for ( unsigned int i = 0; i < nlines; i ++ )
            {
                const kvs::ValueArray<float>& line_coords = line->coords();
                const kvs::ValueArray<unsigned char>& line_colors = line->colors();
                const unsigned int start = line->connection( i ).x();
                const unsigned int end = line->connection( i ).y();
                const unsigned int current = start;

                for ( unsigned int j = 0; j < 3; j ++ )
                {
                    coords[ 3*i+j ] = line_coords[ 3*current+j ];
                    colors[ 3*i+j ] = line_colors[ 3*current+j ];
                }
            }

            point->setCoords( coords );
            point->setColors( colors );
            point->setSize( 1.0f );
            point->setMinMaxObjectCoords( volume->minObjectCoord(), volume->maxObjectCoord() );
            point->setMinMaxExternalCoords( volume->minExternalCoord(), volume->maxExternalCoord() );

            time_step += speed;

            break;
        }
        case kvs::LineObject::Uniline:
        {
            break;
        }
        default:
        {
            std::cout << "Line type not support by filter yet" << std::endl;
            return 1;
        }
    }

    //renderer
    const size_t repeat_level = arg.hasOption( "r" ) ? arg.optionValue<size_t>( "r" ) : 1;
    renderer = new kvs::glew::StochasticRenderer( repeat_level );
    renderer->enableLODControl();

    kvs::NullObject* null = new::kvs::NullObject( volume );

    kvs::TransferFunction tfunc;
    tfunc.setColorMap( kvs::RGBFormulae::Hot( 256 ) );
    kvs::glew::StochasticVolumeRenderer* volume_renderer = new kvs::glew::StochasticVolumeRenderer( volume );
    volume_renderer->setShader( kvs::Shader::BlinnPhong() );
    volume_renderer->setTransferFunction( tfunc );
    volume_renderer->setEdgeSize( 3 );
    volume_renderer->disableShading();
    renderer->registerRenderer( volume_renderer );

    polygon->setOpacity( 20 );
    polygon->setColor( kvs::RGBColor( 255, 255, 255 ) );
    kvs::glew::StochasticPolygonRenderer* polygon_renderer = new kvs::glew::StochasticPolygonRenderer( polygon );
    polygon_renderer->setShader( kvs::Shader::BlinnPhong() );
    renderer->registerRenderer( polygon_renderer );

    point_renderer = new kvs::glew::StochasticPointRenderer( point );
    point_renderer->setPointSize( 1.0f );
    renderer->registerRenderer( point_renderer );

    screen.registerObject( null, renderer );

    return ( app.run() );
}
