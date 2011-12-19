#include "MainScreen.h"
#include "ControlScreen.h"
#include <kvs/glut/TransferFunctionEditor>

class TransferFunctionEditor : public kvs::glut::TransferFunctionEditor
{
private:

    kvs::glew::StochasticRenderer*          m_ref_renderer;
    kvs::glew::StochasticVolumeRenderer*    m_ref_volume_renderer;

public:

    TransferFunctionEditor( 
            kvs::glut::Screen* screen,
            kvs::glew::StochasticRenderer* renderer,
            kvs::glew::StochasticVolumeRenderer* volume_renderer):
            kvs::glut::TransferFunctionEditor( screen ),
            m_ref_renderer( renderer ),
            m_ref_volume_renderer( volume_renderer )
    {
    }

    void apply( void )
    {        
        m_ref_volume_renderer->setTransferFunction( transferFunction() );
        m_ref_renderer->clearEnsembleBuffer();
        screen()->redraw();
    }
};

class Argument : public kvs::CommandLine
{

public:

    Argument( int argc, char** argv ) :
        CommandLine( argc, argv )
    {
        add_help_option();
        add_option( "r", "[size_t] Repeat level. ( default : 1 )", 1, false );

        add_option( "volume", "[string] kvs::UnstructuredVolumeObject file path. ( optional )", 1, false );
        add_option( "polygon", "[string] kvs::PolygonObject file path. ( optional )", 1, false );
        add_option( "tfunc", "[string] kvs::TransferFunction file path. ( optional )", 1, false );

        add_option( "DisableShading", "Disable shading. ( default : eable shading )", 0, false );

        if( !this->parse() ) exit( EXIT_FAILURE );
    }

};

class TimerEvent : public kvs::TimerEventListener
{
public:
    void update ( kvs::TimeEvent* event )
    {
        screen()->redraw();
    }
};

int main( int argc, char** argv )
{
    kvs::glut::Application app( argc, argv );
    Argument arg( app.argc(), app.argv() );

    //volume
    std::string default_volume( "../../data/engine/v6engine.kvsml" );
    kvs::VolumeObjectBase* volume = arg.hasOption( "volume2" ) ?
        new kvs::UnstructuredVolumeImporter( arg.optionValue<std::string>( "volume2" ) ):
        new kvs::UnstructuredVolumeImporter( default_volume );

    if ( !volume )
    {
        kvsMessageError( "Cannot create an unstructured volume object." );
        exit( EXIT_FAILURE );
    }

    kvs::MainScreen main_screen( &app );
    int interval = 30;
    kvs::glut::Timer timer( interval );
    TimerEvent timer_event;
    main_screen.addTimerEvent( &timer_event, &timer );
    main_screen.show();

    TransferFunctionEditor editor( &main_screen, main_screen.renderer(), main_screen.volumeRenderer() );
    editor.setVolumeObject( volume );
    editor.show();

    kvs::ControlScreen control_screen( &app );
    control_screen.attachMainScreen( &main_screen );
    control_screen.show();

    return( app.run() );
}
