#include <kvs/Bounds>
#include <kvs/CommandLine>
#include <kvs/ExternalFaces>
#include <kvs/LineExporter>
#include <kvs/PointRenderer>
#include <kvs/PolygonImporter>
#include <kvs/MouseButton>
#include <kvs/RGBFormulae>
#include <kvs/UnstructuredVectorToScalar>
#include <kvs/UnstructuredVolumeImporter>

#include <kvs/glut/Application>
#include <kvs/glut/CheckBox>
#include <kvs/glut/LegendBar>
#include <kvs/glut/PushButton>
#include <kvs/glut/RadioButtonGroup>
#include <kvs/glut/Screen>
#include <kvs/glut/Slider>
#include <kvs/glut/TransferFunctionEditor>

#include "CellLocatorBIH.h"
//#include "ControlScreen.h"
#include "CubicPointObject.h"
#include "HyperStreamline.h"

#include "PolygonToPolygon.h"

#define USE_KVS

#ifdef USE_KVS
#include <kvs/glew/StochasticRenderingCompositor>
#include <kvs/glew/StochasticLineEngine>
#include <kvs/glew/StochasticPointEngine>
#include <kvs/glew/StochasticPolygonEngine>
#include <kvs/glew/StochasticTetrahedraEngine>
#else
#include "NullObject.h"
#include "PolygonToPolygon.h"
#include "StochasticLineRenderer.h"
#include "StochasticPointRenderer.h"
#include "StochasticPolygonRenderer.h"
#include "StochasticRenderer.h"
#include "StochasticVolumeRenderer.h"
#endif

#ifdef USE_KVS
kvs::glew::StochasticPointEngine*         m_point_renderer        = NULL;
kvs::glew::StochasticLineEngine*          m_line_renderer         = NULL;
kvs::glew::StochasticPolygonEngine*       m_polygon_renderer      = NULL;
kvs::glew::StochasticTetrahedraEngine*    m_volume_renderer       = NULL;
kvs::glew::StochasticRenderingCompositor* m_compositor            = NULL;
#else
kvs::NullObject*                        null                      = NULL;
kvs::glew::StochasticRenderer*          m_renderer                = NULL;
kvs::glew::StochasticPointRenderer*     m_point_renderer          = NULL;
kvs::glew::StochasticLineRenderer*      m_line_renderer           = NULL;
kvs::glew::StochasticPolygonRenderer*   m_polygon_renderer        = NULL;
kvs::glew::StochasticVolumeRenderer*    m_volume_renderer         = NULL;
#endif

kvs::glut::Screen*                      p_main_screen             = NULL;
kvs::UnstructuredVolumeObject*          m_volume1                 = NULL;
kvs::UnstructuredVolumeObject*          m_volume2                 = NULL;
kvs::CellLocatorBIH*                    m_locator                 = NULL;
kvs::CubicPointObject*                  m_seed_point              = NULL;
kvs::PolygonObject*                     m_polygon                 = NULL;
kvs::HyperStreamline*                   m_streamline              = NULL;
kvs::HyperStreamline*                   m_buffered_streamline     = NULL;

kvs::TransferFunction                   m_tfunc;  
std::vector<unsigned int>               m_buffer;  
kvs::Vector3f                           center;

unsigned int opacity_polygon = 30;
unsigned int opacity_line    = 60;
unsigned int opacity_point   = 40;
unsigned int tfunc_counter   = 0;
unsigned int line_direction  = 0;

unsigned int resx;
unsigned int resy;
unsigned int resz;
float xmin;
float xmax;
float ymin;
float ymax;
float zmin;
float zmax;

bool flag_locator_built = false;
bool flag_clear_buffer = true;
bool flag_show_volume = false;
bool flag_show_polygon = true;
bool flag_enable_cache = true;

class Argument : public kvs::CommandLine
{

public:

    Argument( int argc, char** argv ) :
        CommandLine( argc, argv )
    {
        add_help_option();
        add_option( "r", "[size_t] Repeat level. ( default : 1 )", 1, false );

        add_option( "volume1",  "[string] Tensor Volume.", 1, false );
        add_option( "volume2",  "[string] Display Volume.", 1, false );
        add_option( "polygon",  "[string] Display Polygon", 1, false );
        add_option( "celltree", "[string] Celltree", 1, false ); 

        if( !this->parse() ) exit( EXIT_FAILURE );
    }

};
Argument* p_arg  = NULL;

#ifdef USE_KVS
#else
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
        this->setTransferFunction( m_ref_volume_renderer->transferFunction() );
    }

    void apply( void )
    {        
        m_ref_volume_renderer->setTransferFunction( transferFunction() );
        m_ref_renderer->clearEnsembleBuffer();
        screen()->redraw();
    }
};
#endif

class PB_INFO : public kvs::glut::PushButton 
{
public:
    PB_INFO( kvs::glut::Screen* screen = 0 ):
        kvs::glut::PushButton( screen ){};

    void pressed()
    {
        this->setCaption("You don't need to press me.  :)");
    }
};
PB_INFO* p_pb_info = NULL;

// Update Streamline
void update_streamline()
{
    if ( !flag_locator_built )
    {
        delete m_locator;
        m_locator = new kvs::CellLocatorBIH();
        m_locator->setDataSet( m_volume1 );
        m_locator->setMode( kvs::CellLocator::CACHEHALF );
        m_locator->setParallel();
        m_locator->initializeCell();
        m_locator->build();
        flag_locator_built = true;
    }

    kvs::HyperStreamline* streamline = new kvs::HyperStreamline();
    streamline->setGoWithNthEigenVector( line_direction );
    streamline->setIntegrationDirection( kvs::HyperStreamline::BothDirections );
    streamline->setLocator( m_locator->cellTree(), m_volume1 );
    if ( !flag_enable_cache )
        streamline->setDisableCache();
    streamline->setIntegrationMethod( kvs::HyperStreamline::RungeKutta2nd );
    streamline->setTransferFunction( m_tfunc );
    streamline->setSeedPoints( m_seed_point );
    streamline->exec( m_volume1 );
    streamline->calculate_color();

    if ( flag_clear_buffer )
    {
#ifdef USE_KVS
        m_compositor->changeObject( m_streamline, streamline, true );
#else
        m_renderer->changeObject( streamline, m_line_renderer, true );
#endif
        m_streamline = streamline;
    }
    else
    {
            //const kvs::ValueArray<kvs::UInt8>& old_colors = m_streamline->colors();
            const std::vector<float>& old_eigenvalues = m_streamline->eigenValues();
            const kvs::ValueArray<kvs::UInt32>& old_connections = m_streamline->connections();
            const kvs::ValueArray<kvs::Real32>& old_coords = m_streamline->coords();

            //const kvs::ValueArray<kvs::UInt8>& new_colors = streamline->colors();
            const std::vector<float>& new_eigenvalues = streamline->eigenValues();
            const kvs::ValueArray<kvs::UInt32>& new_connections = streamline->connections();
            const kvs::ValueArray<kvs::Real32>& new_coords = streamline->coords();

            //const unsigned int size_color = old_colors.size() + new_colors.size();
            const unsigned int size_eigen = old_eigenvalues.size() + new_eigenvalues.size();
            const unsigned int size_conne = old_connections.size() + new_connections.size();
            const unsigned int offs_conne = old_connections.size() > 1 ? old_connections.back() + 1 : 0;
            const unsigned int size_coord = old_coords.size() + new_coords.size();

            //kvs::ValueArray<kvs::UInt8> colors( size_color );
            std::vector<float> eigenvalues;
            eigenvalues.reserve( size_eigen );
            kvs::ValueArray<kvs::UInt32> connections( size_conne );
            kvs::ValueArray<kvs::Real32> coords( size_coord );

            eigenvalues.insert( eigenvalues.end(), old_eigenvalues.begin(), old_eigenvalues.end() );
            eigenvalues.insert( eigenvalues.end(), new_eigenvalues.begin(), new_eigenvalues.end() );

            //for ( unsigned int i = 0; i < old_colors.size(); i ++ )
            //    colors[ i ] = old_colors[ i ];
            //for ( unsigned int i = 0; i < new_colors.size(); i ++ )
            //    colors[ i + old_colors.size() ] = new_colors[ i ];

            for ( unsigned int i = 0; i < old_connections.size(); i ++ )
                connections[ i ] = old_connections[i];
            for ( unsigned int i = 0; i < new_connections.size(); i ++ )
                connections[ i + old_connections.size() ] = new_connections[i] + offs_conne;

            for ( unsigned int i = 0; i < old_coords.size(); i ++ )
                coords[ i ] = old_coords[ i ];
            for ( unsigned int i = 0; i < new_coords.size(); i ++ )
                coords[ i + old_coords.size() ] = new_coords[i];

            delete m_buffered_streamline;
            m_buffered_streamline = new kvs::HyperStreamline();

            //m_buffered_streamline->setColors( colors );
            m_buffered_streamline->setEigenValues( eigenvalues );
            m_buffered_streamline->calculate_color();
            m_buffered_streamline->setEigenValues( eigenvalues );
            m_buffered_streamline->setConnections( connections );
            m_buffered_streamline->setCoords( coords );
            m_buffered_streamline->setLineType( kvs::LineObject::Polyline );
            m_buffered_streamline->setColorType( kvs::LineObject::VertexColor );

#ifdef USE_KVS
            m_compositor->changeObject( m_streamline, m_buffered_streamline, false );
#else
            m_renderer->changeObject( m_buffered_streamline, m_line_renderer, false );
#endif
            m_streamline = m_buffered_streamline;
    }


}
class PB04 : public kvs::glut::PushButton
{
public:
    PB04 ( kvs::glut::Screen* screen = 0 ) :
        kvs::glut::PushButton( screen ){}

    void pressed()
    {
        update_streamline();
    }

};
PB04* p_pb04 = NULL;

// dx
class SLD0 : public kvs::glut::Slider
{
public:

    SLD0( kvs::glut::Screen* screen = 0 ):
        kvs::glut::Slider( screen ){};

    void sliderMoved( void )
    {
        const float dx_new = this->value();
        const float _xmin = center.x() - dx_new / 2;
        const float _xmax = center.x() + dx_new / 2;

        m_seed_point->setXMin( _xmin );
        m_seed_point->setXMax( _xmax );
        m_seed_point->reset_coordinates();

#ifdef USE_KVS
        m_compositor->changeObject( m_seed_point, m_seed_point, false );
#else
        m_renderer->changeObject( m_seed_point, m_point_renderer, false );
#endif
    }
};
SLD0* p_sld0 = NULL;

// dy
class SLD1 : public kvs::glut::Slider
{
public:

    SLD1( kvs::glut::Screen* screen = 0 ):
        kvs::glut::Slider( screen ){};

    void sliderMoved( void )
    {
        const float dy_new = this->value();
        const float _ymin = center.y() - dy_new / 2;
        const float _ymax = center.y() + dy_new / 2;

        m_seed_point->setYMin( _ymin );
        m_seed_point->setYMax( _ymax );
        m_seed_point->reset_coordinates();

#ifdef USE_KVS
        m_compositor->changeObject( m_seed_point, m_seed_point, false );
#else
        m_renderer->changeObject( m_seed_point, m_point_renderer, false );
#endif
    }
};
SLD1* p_sld1 = NULL;

// dz
class SLD2 : public kvs::glut::Slider
{
public:

    SLD2( kvs::glut::Screen* screen = 0 ):
        kvs::glut::Slider( screen ){}

    void sliderMoved( void )
    {
        const float dz_new = this->value();
        const float _zmin = center.z() - dz_new / 2;
        const float _zmax = center.z() + dz_new / 2;

        m_seed_point->setZMin( _zmin );
        m_seed_point->setZMax( _zmax );
        m_seed_point->reset_coordinates();

#ifdef USE_KVS
        m_compositor->changeObject( m_seed_point, m_seed_point, false );
#else
        m_renderer->changeObject( m_seed_point, m_point_renderer, false );
#endif
    }
};
SLD2* p_sld2 = NULL;

// resx 
class SLD3 : public kvs::glut::Slider
{
public:

    SLD3( kvs::glut::Screen* screen = 0 ):
        kvs::glut::Slider( screen ){}

    void sliderMoved( void )
    {
        const unsigned int resx_new = static_cast<int>( this->value() + 0.5f );

        m_seed_point->setResX( resx_new );
        m_seed_point->reset_coordinates();

#ifdef USE_KVS
        m_compositor->changeObject( m_seed_point, m_seed_point, false );
#else
        m_renderer->changeObject( m_seed_point, m_point_renderer, false );
#endif
    }
};
SLD3* p_sld3 = NULL;

// resy 
class SLD4 : public kvs::glut::Slider
{
public:

    SLD4( kvs::glut::Screen* screen = 0 ):
        kvs::glut::Slider( screen ){};

    void sliderMoved( void )
    {
        const unsigned int resy_new = static_cast<int>( this->value() + 0.5f );

        m_seed_point->setResY( resy_new );
        m_seed_point->reset_coordinates();

#ifdef USE_KVS
        m_compositor->changeObject( m_seed_point, m_seed_point, false );
#else
        m_renderer->changeObject( m_seed_point, m_point_renderer, false );
#endif
    }
};
SLD4* p_sld4 = NULL;

// resz 
class SLD5 : public kvs::glut::Slider
{
public:

    SLD5( kvs::glut::Screen* screen = 0 ):
        kvs::glut::Slider( screen ){};

    void sliderMoved( void )
    {
        const unsigned int resz_new = static_cast<int>( this->value() + 0.5f );

        m_seed_point->setResZ( resz_new );
        m_seed_point->reset_coordinates();

#ifdef USE_KVS
        m_compositor->changeObject( m_seed_point, m_seed_point, false );
#else
        m_renderer->changeObject( m_seed_point, m_point_renderer, false );
#endif
    }
};
SLD5* p_sld5 = NULL;

// locator (Optional)
class PB03 : public kvs::glut::PushButton 
{
public:
    PB03 ( kvs::glut::Screen* screen = 0 ) :
        kvs::glut::PushButton( screen ){}

    void pressed()
    {
        std::string default_locator( "celltree.dat" );
        if ( p_arg->hasOption( "celltree" ) )
            m_locator->read( p_arg->optionValue<std::string>( "celltree" ) );
        else
            m_locator->read( default_locator );
        
        flag_locator_built = true;
    }
};
PB03* p_pb03 = NULL;

// polygon (External)
class PB02 : public kvs::glut::PushButton 
{
public:
    PB02 ( kvs::glut::Screen* screen = 0 ) :
        kvs::glut::PushButton( screen ){}

    void pressed()
    {
        flag_show_polygon = !flag_show_polygon;
        if ( !flag_show_polygon )
        {
            //empty polygon 
            std::vector<kvs::Real32> coords;
            std::vector<kvs::UInt32> connections;
            std::vector<float>       values;
            coords.push_back( 0.0 );
            coords.push_back( 0.0 );
            coords.push_back( 0.0 );
            for ( unsigned int i = 0; i < 3; i ++ )
                connections.push_back( 0 );
            values.push_back(0);

            kvs::PolygonObject* temp_polygon = new kvs::PolygonObject();
            temp_polygon->setCoords( kvs::ValueArray<kvs::Real32>(coords) );
            temp_polygon->setConnections( kvs::ValueArray<kvs::UInt32>(connections) );
            temp_polygon->setColor( kvs::RGBColor( 255,255,255 ) );
            temp_polygon->setOpacity( 0 );
            temp_polygon->setPolygonType( kvs::PolygonObject::Triangle );

#ifdef USE_KVS
            m_compositor->changeObject( m_polygon, temp_polygon, true );
#else
            m_renderer->changeObject( temp_polygon, m_polygon_renderer, true );
#endif
            m_polygon = temp_polygon;

            this->setCaption( "Show Polygon" );
            p_pb_info->setCaption( "Polygon Read!" );
        }
        else
        {
#ifdef WIN32
            std::string default_polygon( "D:\\Koron\\Dropbox\\Work\\Viz\\Hyperstreamline\\data\\engine\\v6engine_external_face.kvsml" );
#else
            std::string default_polygon( "../../data/engine/v6engine_external_face.kvsml" );
#endif
            kvs::PolygonObject* temp_polygon = NULL;
            if ( p_arg->hasOption( "polygon" ) )
               temp_polygon = new kvs::PolygonImporter( p_arg->optionValue<std::string>("polygon") );
            else
            {
                kvs::PolygonObject* import_polygon = new kvs::PolygonImporter( default_polygon );
                temp_polygon = new kvs::PolygonToPolygon( import_polygon );
            }

            if ( !temp_polygon )
            {
                kvsMessageError( "Cannot create surface." );
                exit( EXIT_FAILURE );
            }

            temp_polygon->setOpacity( opacity_polygon );
            //temp_polygon->setColor( kvs::RGBColor( 255,255,255 ) );

#ifdef USE_KVS
            m_compositor->changeObject( m_polygon, temp_polygon, true );
#else
            m_renderer->changeObject( temp_polygon, m_polygon_renderer, true );
#endif
            m_polygon = temp_polygon;

            this->setCaption( "Hide Polygon" );
            p_pb_info->setCaption( "Polygon Hided!" );
        }
        this->deactivate();
    }

};
PB02* p_pb02 = NULL;

// volume2 (Displacement)
class PB01 : public kvs::glut::PushButton 
{
public:
    PB01 ( kvs::glut::Screen* screen = 0 ) :
        kvs::glut::PushButton( screen ){}

    void pressed()
    {
        flag_show_volume = !flag_show_volume;
        if ( flag_show_volume )
        {
#ifdef WIN32
            std::string default_volume2( "D:\\Koron\\Dropbox\\Work\\Viz\\Hyperstreamline\\data\\engine\\v6engine_displacement.kvsml" );
#else
            std::string default_volume2( "../../data/engine/v6engine_displacement.kvsml" );
#endif
            kvs::UnstructuredVolumeObject* temp_volume2 = NULL;
            if ( p_arg->hasOption( "volume2" ) )
                temp_volume2 = new kvs::UnstructuredVolumeImporter( p_arg->optionValue<std::string>("volume2") );
            else
                temp_volume2 = new kvs::UnstructuredVolumeImporter( default_volume2 );
                
            if ( temp_volume2->veclen() == 3 )
            {
                kvs::UnstructuredVolumeObject* volume2 = new kvs::UnstructuredVectorToScalar( temp_volume2 );
                delete temp_volume2;
                temp_volume2 = volume2;
            }  

#ifdef USE_KVS
            m_compositor->changeObject( m_volume2, temp_volume2, true );
#else
            m_renderer->changeObject( temp_volume2, m_volume_renderer, true );
#endif
            m_volume2 = temp_volume2;

            p_pb_info->setCaption( "Displacement Volume Read!" );
            this->setCaption( "Hide Volume" );
        }
        else
        {
            //empty volume
            std::vector<kvs::Real32> coords;
            std::vector<kvs::UInt32> connections;
            std::vector<float>       values;
            coords.push_back( 0.0 );
            coords.push_back( 0.0 );
            coords.push_back( 0.0 );
            for ( unsigned int i = 0; i < 4; i ++ )
                connections.push_back( 0 );
            values.push_back(0);

            kvs::UnstructuredVolumeObject* volume2 = new kvs::UnstructuredVolumeObject();
            volume2->setCellType( kvs::UnstructuredVolumeObject::Tetrahedra );
            volume2->setConnections( kvs::ValueArray<kvs::UInt32>( connections ) );
            volume2->setCoords( kvs::ValueArray<float>( coords ) );
            volume2->setValues( kvs::AnyValueArray( values ) );
            volume2->setNCells(0);
            volume2->setNNodes(0);

#ifdef USE_KVS
            m_compositor->changeObject( m_volume2, volume2, true );
#else
            m_renderer->changeObject( volume2, m_volume_renderer, true );
#endif
            m_volume2 = volume2;

            p_pb_info->setCaption( "Volume Hided!" );
            this->setCaption( "Show Volume" );
        }
     }
};
PB01* p_pb01 = NULL;

// volume1 (Tensor)
class PB00 : public kvs::glut::PushButton 
{
public:
    PB00 ( kvs::glut::Screen* screen = 0 ) :
        kvs::glut::PushButton( screen ){}

    void pressed()
    {
#ifdef WIN32
        std::string default_volume1( "D:\\Koron\\Dropbox\\Work\\Viz\\Hyperstreamline\\data\\engine\\v6engine_stress_and_mises.kvsml" );
#else
        std::string default_volume1( "../../data/engine/v6engine_stress_and_mises.kvsml" );
#endif
        if ( p_arg->hasOption( "volume1" ) )
            m_volume1 = new kvs::UnstructuredVolumeImporter( p_arg->optionValue<std::string>("volume1") );
        else
            m_volume1 = new kvs::UnstructuredVolumeImporter( default_volume1 );

        if ( m_volume1->veclen() != 6 && m_volume1->veclen() != 7 )
        {
            p_pb_info->setCaption( "Veclen wrong!" );
            return;
        }

        // initialize sliders
        const kvs::Vector3f min_coord = m_volume1->minExternalCoord();
        const kvs::Vector3f max_coord = m_volume1->maxExternalCoord();
        xmin = min_coord.x();
        ymin = min_coord.y();
        zmin = min_coord.z();
        xmax = max_coord.x();
        ymax = max_coord.y();
        zmax = max_coord.z();
        center.x() = ( xmin + xmax ) / 2;
        center.y() = ( ymin + ymax ) / 2;
        center.z() = ( zmin + zmax ) / 2;
        resx = 10;
        resy = 10;
        resz = 10;
        const float dx = xmax - xmin;
        const float dy = ymax - ymin;
        const float dz = zmax - zmin;

        p_sld0->setValue( dx );
        p_sld0->setRange( 0.0, dx );
        p_sld0->activate();
        p_sld1->setValue( dy );
        p_sld1->setRange( 0.0, dy );
        p_sld1->activate();
        p_sld2->setValue( dz );
        p_sld2->setRange( 0, dz );
        p_sld2->activate();

        p_sld3->setValue( resx );
        p_sld3->setRange( 0, 30 );
        p_sld3->activate();
        p_sld4->setValue( resy );
        p_sld4->setRange( 0, 30 );
        p_sld4->activate();
        p_sld5->setValue( resz );
        p_sld5->setRange( 0, 30 );
        p_sld5->activate();

        m_seed_point->setResX( resx );
        m_seed_point->setResY( resy );
        m_seed_point->setResZ( resz );
        m_seed_point->setXMin( xmin );
        m_seed_point->setXMax( xmax );
        m_seed_point->setYMin( ymin );
        m_seed_point->setYMax( ymax );
        m_seed_point->setZMin( zmin );
        m_seed_point->setZMax( zmax );
        m_seed_point->reset_coordinates();
        m_seed_point->setMinMaxObjectCoords( m_volume1->minObjectCoord(), m_volume1->maxObjectCoord() );
        m_seed_point->setMinMaxExternalCoords( m_volume1->minExternalCoord(), m_volume1->maxExternalCoord() );
        
#ifdef USE_KVS
#else
        m_renderer->changeObject( m_seed_point, m_point_renderer, false );
#endif

        p_pb_info->setCaption( "Tensor Volume Read!" );
        std::cout << "Tensor Volume Read" << std::endl;
        this->deactivate();
    }

};
PB00* p_pb00 = NULL;

// enable line buffer
class PB05 : public kvs::glut::PushButton 
{
public:
    PB05 ( kvs::glut::Screen* screen = 0 ) :
        kvs::glut::PushButton( screen ){}

    void pressed()
    {
        flag_clear_buffer = false;
        p_pb_info->setCaption("Line Buffer Mode enabled");
    }
};
PB05* p_pb05 = NULL;

// clear
class PB06 : public kvs::glut::PushButton 
{
public:
    PB06 ( kvs::glut::Screen* screen = 0 ) :
        kvs::glut::PushButton( screen ){}

    void pressed()
    {
        flag_clear_buffer = true;

        std::vector<kvs::Real32> coords;
        std::vector<kvs::UInt32> connections;
        std::vector<kvs::UInt8>  colors;
        std::vector<kvs::Real32>  eigenvalue;
        coords.push_back( 0.0 );
        coords.push_back( 0.0 );
        coords.push_back( 0.0 );
        connections.push_back( 0 );
        connections.push_back( 0 );
        colors.push_back( 0 );
        colors.push_back( 0 );
        colors.push_back( 0 );
        eigenvalue.push_back(0);

        kvs::HyperStreamline* streamline = new kvs::HyperStreamline();
        streamline->setLineType( kvs::LineObject::Polyline );
        streamline->setColorType( kvs::LineObject::VertexColor );
        streamline->setCoords( kvs::ValueArray<kvs::Real32>( coords ) );
        streamline->setConnections( kvs::ValueArray<kvs::UInt32>( connections ) );
        streamline->setColors( kvs::ValueArray<kvs::UInt8>( colors ) );
        streamline->setEigenValues( eigenvalue );
        streamline->setSize( 1.0f );

        delete m_buffered_streamline;
        m_buffered_streamline = NULL;

#ifdef USE_KVS
        m_compositor->changeObject( m_streamline, streamline, false );
#else
        m_renderer->changeObject( streamline, m_line_renderer, false );
#endif
        m_streamline = streamline;

        p_pb_info->setCaption("All in a mess huh? let's do it again");

    }
};
PB06* p_pb06 = NULL;

// save celltree
class PB07 : public kvs::glut::PushButton 
{
public:
    PB07 ( kvs::glut::Screen* screen = 0 ) :
        kvs::glut::PushButton( screen ){}

    void pressed()
    {
        // locator (Optional)
        if ( m_locator == NULL )
        {
            p_pb_info->setCaption("CellTree hasn't been built yet!!");
            return;
        }

        std::string default_locator( "celltree.dat" );
        m_locator->write( default_locator );

        p_pb_info->setCaption("CellTree saved as celltree.dat");
    }
};
PB07* p_pb07 = NULL;

// major principal stress direction
class RB00 : public kvs::glut::RadioButton
{
public:
    RB00( kvs::glut::Screen* screen ):
      kvs::glut::RadioButton( screen ){}

    void stateChanged()
    {
        if ( this->state() )
        {
            line_direction = 0;
        }
    }
};

// intermediate principal stress direction
class RB01 : public kvs::glut::RadioButton
{
public:
    RB01( kvs::glut::Screen* screen ):
      kvs::glut::RadioButton( screen ){}

    void stateChanged()
    {
        if ( this->state() )
        {
            line_direction = 1;
        }
    }
};

// minor principal stress direction
class RB02 : public kvs::glut::RadioButton
{
public:
    RB02( kvs::glut::Screen* screen ):
      kvs::glut::RadioButton( screen ){}

    void stateChanged()
    {
        if ( this->state() )
        {
            line_direction = 2;
        }
    }
};

class CB00 : public kvs::glut::CheckBox
{
public:

    CB00( kvs::glut::Screen* screen ):
        kvs::glut::CheckBox( screen ){};

    void stateChanged( void )
    {
        flag_enable_cache = this->state();
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

class PaintEvent : public kvs::PaintEventListener
{
public:

    PaintEvent() : kvs::PaintEventListener(){}

    void update()
    {
        // draw the controllable area
         GLint vp[4]; glGetIntegerv( GL_VIEWPORT, vp );
         const GLint left   = vp[0];
         const GLint bottom = vp[1];
         const GLint right  = vp[2];
         const GLint top    = vp[3];
     
         glPushAttrib( GL_ALL_ATTRIB_BITS );
         glMatrixMode( GL_MODELVIEW );  glPushMatrix(); glLoadIdentity();
         glMatrixMode( GL_PROJECTION ); glPushMatrix(); glLoadIdentity();
         glOrtho( left, right, top, bottom, -1, 1 ); // The origin is upper-left.
         glDisable( GL_DEPTH_TEST );        

        const GLfloat x0 = 20;
        const GLfloat y0 = 20;
        const GLfloat x1 = 330;
        const GLfloat y1 = 400;
        const GLubyte r = 255;
        const GLubyte g = 0;
        const GLubyte b = 0;

        glColor3ub( r, g, b );
        glLineWidth( 1 );

        glBegin( GL_LINE_STRIP );
            glVertex2f( x0, y1 );
            glVertex2f( x0, y0 );
            glVertex2f( x1, y0 );
            glVertex2f( x1, y1 );
            glVertex2f( x0, y1 );
        glEnd();    
            
        glPopMatrix();
        glMatrixMode( GL_MODELVIEW );
        glPopMatrix();
        glPopAttrib();
    }
};

class KeyPressEvent : public kvs::KeyPressEventListener
{
public:

    void update( kvs::KeyEvent* event )
    {
        switch ( event->key() )
        {
            case kvs::Key::o: this->screen()->controlTarget() = kvs::ScreenBase::TargetObject; break;
            case kvs::Key::l: this->screen()->controlTarget() = kvs::ScreenBase::TargetLight; break;
            case kvs::Key::Up:
            {
                if ( opacity_polygon < 255 )
                    opacity_polygon += 5;
                m_polygon->setOpacity( opacity_polygon );
#ifdef USE_KVS
                m_compositor->clearEnsembleBuffer();
#else
                m_renderer->clearEnsembleBuffer();
#endif
                this->screen()->redraw();
                break;
            }
            case kvs::Key::Down:
            {
                if ( opacity_polygon > 0 )
                    opacity_polygon -= 5;
                m_polygon->setOpacity( opacity_polygon );
#ifdef USE_KVS
                m_compositor->clearEnsembleBuffer();
#else
                m_renderer->clearEnsembleBuffer();
#endif
                this->screen()->redraw();
                break;
            }
            case kvs::Key::Left:
            {
                if ( opacity_line < 255 )
                    opacity_line += 5;
                m_line_renderer->setOpacity( opacity_line );
#ifdef USE_KVS
                m_compositor->clearEnsembleBuffer();
#else
                m_renderer->clearEnsembleBuffer();
#endif
                this->screen()->redraw();
                break;
            }
            case kvs::Key::Right:
            {
                if ( opacity_line > 0 )
                    opacity_line -= 5;
                m_line_renderer->setOpacity( opacity_line );
#ifdef USE_KVS
                m_compositor->clearEnsembleBuffer();
#else
                m_renderer->clearEnsembleBuffer();
#endif
                this->screen()->redraw();
                break;
            }
#ifdef UNIFORM_DISTRIBUTION
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
#endif
            case kvs::Key::t:
            {
                switch( tfunc_counter )
                {
                    case 0: m_tfunc.setColorMap( kvs::RGBFormulae::AFMHot( 256 ) ); break;
                    case 1: m_tfunc.setColorMap( kvs::RGBFormulae::Bone( 256 ) ); break;
                    case 2: m_tfunc.setColorMap( kvs::RGBFormulae::Hot( 256 ) ); break;
                    case 3: m_tfunc.setColorMap( kvs::RGBFormulae::Jet( 256 ) ); break;
                    case 4: m_tfunc.setColorMap( kvs::RGBFormulae::Ocean( 256 ) ); break;
                    case 5: m_tfunc.setColorMap( kvs::RGBFormulae::PM3D( 256 ) ); break;
                    case 6: m_tfunc.setColorMap( kvs::RGBFormulae::Rainbow( 256 ) ); break;
                    default : break;
                }

                tfunc_counter++;
                if ( tfunc_counter > 6 )
                    tfunc_counter = 0;

                m_volume_renderer->setTransferFunction( m_tfunc );
#ifdef USE_KVS
                m_compositor->clearEnsembleBuffer();
#else
                m_renderer->clearEnsembleBuffer();
#endif
                this->screen()->redraw();
                break;
            }
            default: break;
        }

    }
};

class MouseMoveEvent : public kvs::MouseMoveEventListener
{
public:

    MouseMoveEvent() : kvs::MouseMoveEventListener()
    {
    }

    void update ( kvs::MouseEvent* event = 0 )
    {
        if ( event->x() > 20 && event->x() < 330 && event->y() > 20 && event->y() < 400 && m_seed_point->nvertices() > 0 )
        {
            const kvs::Xform x = p_main_screen->objectManager()->xform();

            const float* pcoord = m_seed_point->coords().pointer();
            const unsigned int nvertices = m_seed_point->nvertices();
            kvs::ValueArray<float> coords( nvertices * 3 );

            if ( event->button() == kvs::MouseButton::Right )
            {
                this->screen()->mouse()->setMode( kvs::Mouse::Translation );
                this->screen()->mouse()->move( event->x(), event->y() );
                kvs::Vector3f translation = this->screen()->mouse()->translation();
                const kvs::Vector3f normalize = p_main_screen->objectManager()->normalize();

                translation.x() /= normalize.x() * x.scaling().x();
                translation.y() /= normalize.y() * x.scaling().y();
                translation.z() /= normalize.z() * x.scaling().z();

                for ( unsigned int i = 0; i < nvertices; i ++ )
                {
                    kvs::Vector3f coord( pcoord );
                    const kvs::Vector3f new_coord = coord + translation * x.rotation();
                    coords[ 3 * i] = new_coord.x();
                    coords[ 3 * i + 1] = new_coord.y();
                    coords[ 3 * i + 2] = new_coord.z();
                    pcoord += 3;
                }
                m_seed_point->setCoords( coords );
                
            }

            if ( event->button() == kvs::MouseButton::Left )
            {
                this->screen()->mouse()->setMode( kvs::Mouse::Rotation );
                this->screen()->mouse()->move( event->x(), event->y() );
                kvs::Matrix33f rotation = this->screen()->mouse()->rotation().toMatrix();

                for ( unsigned int i = 0; i < nvertices; i ++ )
                {
                    kvs::Vector3f coord( pcoord );
                    const kvs::Vector3f new_coord = coord * rotation;
                    coords[ 3 * i] = new_coord.x();
                    coords[ 3 * i + 1] = new_coord.y();
                    coords[ 3 * i + 2] = new_coord.z();
                    pcoord += 3;
                }
                m_seed_point->setCoords( coords );
            }
#ifdef USE_KVS
            m_compositor->changeObject( m_seed_point, m_seed_point, false );
#else
            m_renderer->changeObject( m_seed_point, m_point_renderer, false );
#endif
            m_seed_point->updateMinMaxCoords();
            const kvs::Vector3f min = m_seed_point->minObjectCoord();
            const kvs::Vector3f max = m_seed_point->maxObjectCoord();
            center = ( min + max ) / 2;
            m_seed_point->setXMin( min.x() );
            m_seed_point->setXMax( max.x() );
            m_seed_point->setYMin( min.y() );
            m_seed_point->setYMax( max.y() );
            m_seed_point->setZMin( min.z() );
            m_seed_point->setZMax( max.z() );

            if ( m_seed_point->nvertices() <= 8 )
                update_streamline();
        }
    }
};

int main( int argc, char** argv )
{
    kvs::glut::Application app( argc, argv );
    Argument arg( app.argc(), app.argv() );
    p_arg = &arg;

    // for empty objects
    std::vector<kvs::Real32> coords;
    std::vector<kvs::UInt32> connections;
    std::vector<kvs::UInt8>  colors;
    std::vector<float>       values;
    std::vector<float>       eigenvalues;

    colors.push_back(255);
    colors.push_back(255);
    colors.push_back(255);
    for ( unsigned int i = 0; i < 3; i++ )
    {
        coords.push_back(0.0);
    }
    for ( unsigned int i = 0; i < 2; i ++ )
    {
        connections.push_back(0);
        values.push_back(0.0);
    }
    eigenvalues.push_back(0);

    // Volume 2 (Displacement)
    m_volume2 = new kvs::UnstructuredVolumeObject();
    m_volume2->setCellType( kvs::UnstructuredVolumeObject::Tetrahedra );
    m_volume2->setConnections( kvs::ValueArray<kvs::UInt32>( connections ) );
    m_volume2->setCoords( kvs::ValueArray<float>( coords ) );
    m_volume2->setValues( kvs::AnyValueArray( values ) );
    m_volume2->setNCells(0);
    m_volume2->setNNodes(0);

    // polygon (Empty External)
#ifdef WIN32
    std::string default_polygon( "D:\\Koron\\Dropbox\\Work\\Viz\\Hyperstreamline\\data\\engine\\v6engine_external_face.kvsml" );
#else
    std::string default_polygon( "../../data/engine/v6engine_external_face.kvsml" );
#endif
    if ( p_arg->hasOption( "polygon" ) )
    {
       kvs::PolygonObject* import_polygon = new kvs::PolygonImporter( p_arg->optionValue<std::string>("polygon") );
       m_polygon = new kvs::PolygonToPolygon( import_polygon );
    }
    else
    {
        kvs::PolygonObject* import_polygon = new kvs::PolygonImporter( default_polygon );
        m_polygon = new kvs::PolygonToPolygon( import_polygon );
    }

    if ( !m_polygon )
    {
        kvsMessageError( "Cannot create surface." );
        exit( EXIT_FAILURE );
    }
    m_polygon->setOpacity( opacity_polygon );
    m_polygon->setColor( kvs::RGBColor( 255,255,255 ) );

    const kvs::Vector3f min_coord = m_polygon->minExternalCoord();
    const kvs::Vector3f max_coord = m_polygon->maxExternalCoord();
    xmin = min_coord.x();
    ymin = min_coord.y();
    zmin = min_coord.z();
    xmax = max_coord.x();
    ymax = max_coord.y();
    zmax = max_coord.z();
    resx = 10;
    resy = 10;
    resz = 10;

    m_seed_point = new kvs::CubicPointObject();
    m_seed_point->reset_coordinates( resx, resy, resz, xmin, xmax, ymin, ymax, zmin, zmax );
    m_seed_point->updateMinMaxCoords();

    // streamlines (Empty LineObject)
    m_streamline = new kvs::HyperStreamline();
    m_streamline->setLineType( kvs::LineObject::Polyline );
    m_streamline->setColorType( kvs::LineObject::VertexColor );
    m_streamline->setCoords( kvs::ValueArray<kvs::Real32>( coords ) );
    m_streamline->setConnections( kvs::ValueArray<kvs::UInt32>( connections ) );
    m_streamline->setEigenValues( eigenvalues );
    m_streamline->setColors( kvs::ValueArray<kvs::UInt8>( colors ) );
    m_streamline->setSize( 1.0f );

    // main screen
    kvs::glut::Screen main_screen( &app );
    p_main_screen = &main_screen;
    //main_screen.background()->setColor( kvs::RGBAColor( 0, 0, 0, 1 ) );
    int interval = 30;
    kvs::glut::Timer timer( interval );
    TimerEvent timer_event;
    KeyPressEvent keypress_event;

    main_screen.addTimerEvent( &timer_event, &timer );
    main_screen.addKeyPressEvent( &keypress_event );
    main_screen.show();

    // renderer
#ifdef USE_KVS
    m_compositor = new kvs::glew::StochasticRenderingCompositor( &main_screen );
    m_compositor->enableLODControl();

    m_tfunc.adjustRange( m_volume2 );
    m_tfunc.setColorMap( kvs::RGBFormulae::AFMHot(256) );

    m_volume_renderer = new kvs::glew::StochasticTetrahedraEngine();
    m_volume_renderer->setTransferFunction( m_tfunc );
    m_volume_renderer->setShader( kvs::Shader::BlinnPhong() );
    m_volume_renderer->disableShading();
    m_volume_renderer->setEdgeSize( 2 );

    m_polygon_renderer = new kvs::glew::StochasticPolygonEngine();
    m_polygon_renderer->setShader( kvs::Shader::BlinnPhong() );
   
    m_point_renderer = new kvs::glew::StochasticPointEngine();
    m_point_renderer->disableShading();

    m_line_renderer = new kvs::glew::StochasticLineEngine();
    m_line_renderer->setShader( kvs::Shader::BlinnPhong() );
    m_line_renderer->setOpacity( opacity_line );

    m_compositor->registerObject( m_volume2, m_volume_renderer );
    m_compositor->registerObject( m_polygon, m_polygon_renderer );
    m_compositor->registerObject( m_seed_point, m_point_renderer );
    m_compositor->registerObject( m_streamline, m_line_renderer );
#else
    m_renderer = new kvs::glew::StochasticRenderer( 1 );
    m_renderer->enableLODControl();

    m_tfunc.adjustRange( m_volume2 );
    m_tfunc.setColorMap( kvs::RGBFormulae::AFMHot(256) );
    m_volume_renderer = new kvs::glew::StochasticVolumeRenderer( m_volume2 );
    m_volume_renderer->setTransferFunction( m_tfunc );
    m_volume_renderer->setShader( kvs::Shader::BlinnPhong() );
    m_volume_renderer->disableShading();
    m_volume_renderer->setEdgeSize( 2 );
    m_renderer->registerRenderer( m_volume_renderer );

    m_polygon_renderer = new kvs::glew::StochasticPolygonRenderer( m_polygon );
    m_polygon_renderer->setShader( kvs::Shader::BlinnPhong() );
    m_renderer->registerRenderer( m_polygon_renderer );

    m_point_renderer = new kvs::glew::StochasticPointRenderer( m_seed_point );
    m_point_renderer->disableShading();
    m_renderer->registerRenderer( m_point_renderer );

    m_line_renderer = new kvs::glew::StochasticLineRenderer( m_streamline );
    m_line_renderer->setShader( kvs::Shader::BlinnPhong() );
    m_line_renderer->setOpacity( opacity_line );
    m_renderer->registerRenderer( m_line_renderer );

    null = new::kvs::NullObject( m_seed_point );
    null->setName( "null" );
    main_screen.registerObject( null, m_renderer );
    main_screen.show();

#endif

    // tfunc editor
    //TransferFunctionEditor editor( &main_screen, m_renderer, m_volume_renderer );
    //editor.setVolumeObject( m_volume2 );
    //editor.show();


    //kvs::ControlScreen control_screen( &app );
    MouseMoveEvent mouse_move_event;
    kvs::glut::Screen control_screen( &app );
    //control_screen.addMouseMoveEvent( &mouse_move_event );
    control_screen.setMouseMoveEvent( &mouse_move_event );
    control_screen.setTitle( "kvs::ControlScreen" );
    control_screen.setGeometry( 512, 0, 600, 560 );

#ifdef USE_KVS
    //control_screen.attachMainScreen( p_main_screen, m_seed_point, m_compositor, m_point_renderer );
#else
    //control_screen.attachMainScreen( p_main_screen, m_seed_point, m_renderer, m_point_renderer );
#endif
    control_screen.show();

    PaintEvent paint_event;
    control_screen.addPaintEvent( &paint_event );

    const int width = control_screen.width();
    const int height = control_screen.height();
    const int ui_width = 240;

    PB_INFO pb_info( &control_screen );
    p_pb_info = &pb_info;
    pb_info.setX( width/6 );
    pb_info.setY( height - 40 );
    pb_info.setWidth( 400 );
    pb_info.setMargin( 5 );
    pb_info.setTextMargin( 5 );
    pb_info.setCaption("");
    //pb_info.deactivate();
    pb_info.show();


    PB00 pb00( &control_screen );
    p_pb00 = &pb00;
    pb00.setX( width - ui_width );
    pb00.setY( 10 );
    pb00.setWidth( 220 );
    pb00.setMargin( 5 );
    pb00.setTextMargin( 5 );
    pb00.setCaption("Read volume1(Stress)");
    pb00.show();

    PB01 pb01( &control_screen );
    p_pb01 = &pb01;
    pb01.setX( pb00.x() );
    pb01.setY( pb00.y() + pb00.height() );
    pb01.setWidth( 220 );
    pb01.setMargin( 5 );
    pb01.setTextMargin( 5 );
    pb01.setCaption("Read volume2(Displacement)");
    pb01.deactivate();
    pb01.show();

    PB02 pb02( &control_screen );
    p_pb02 = &pb02;
    pb02.setX( pb01.x() );
    pb02.setY( pb01.y() + pb01.height() );
    pb02.setWidth( 220 );
    pb02.setMargin( 5 );
    pb02.setTextMargin( 5 );
    pb02.setCaption("Hide Polygon");
    //pb02.deactivate();
    pb02.show();

    PB03 pb03( &control_screen );
    p_pb03 = &pb03;
    pb03.setX( pb02.x() );
    pb03.setY( pb02.y() + pb02.height() );
    pb03.setWidth( 220 );
    pb03.setMargin( 5 );
    pb03.setTextMargin( 5 );
    pb03.setCaption("Read CellTree(Optional)");
    //pb03.deactivate();
    pb03.show();

    SLD0 sld0( &control_screen );
    p_sld0 = &sld0;
    sld0.setX( pb03.x() );
    sld0.setY( pb03.y() + pb03.height() );
    sld0.setWidth( 220 );
    sld0.setMargin( 0 );
    sld0.setCaption("dx");
    //sld0.deactivate();
    sld0.show();

    SLD1 sld1( &control_screen );
    p_sld1 = &sld1;
    sld1.setX( sld0.x() );
    sld1.setY( sld0.y() + sld0.height() );
    sld1.setWidth( 220 );
    sld1.setMargin( 0 );
    sld1.setCaption("dy");
    //sld1.deactivate();
    sld1.show();

    SLD2 sld2( &control_screen );
    p_sld2 = &sld2;
    sld2.setX( sld1.x() );
    sld2.setY( sld1.y() + sld1.height() );
    sld2.setWidth( 220 );
    sld2.setMargin( 0 );
    sld2.setCaption("dz");
    //sld2.deactivate();
    sld2.show();

    SLD3 sld3( &control_screen );
    p_sld3 = &sld3;
    sld3.setX( sld2.x() );
    sld3.setY( sld2.y() + sld2.height() );
    sld3.setWidth( 220 );
    sld3.setMargin( 0 );
    sld3.setCaption("resx");
    //sld3.deactivate();
    sld3.show();

    SLD4 sld4( &control_screen );
    p_sld4 = &sld4;
    sld4.setX( sld3.x() );
    sld4.setY( sld3.y() + sld3.height() );
    sld4.setWidth( 220 );
    sld4.setMargin( 0 );
    sld4.setCaption("resy");
    //sld4.deactivate();
    sld4.show();

    SLD5 sld5( &control_screen );
    p_sld5 = &sld5;
    sld5.setX( sld4.x() );
    sld5.setY( sld4.y() + sld4.height() );
    sld5.setWidth( 220 );
    sld5.setMargin( 0 );
    sld5.setCaption("resz");
    //sld5.deactivate();
    sld5.show();

    PB04 pb04( &control_screen );
    p_pb04 = &pb04;
    pb04.setX( sld5.x() );
    pb04.setY( sld5.y() + sld5.height() );
    pb04.setWidth( 220 );
    pb04.setMargin( 15 );
    pb04.setTextMargin( 5 );
    pb04.setCaption("Update Streamline");
    //pb04.deactivate();
    pb04.show();

    PB05 pb05( &control_screen );
    p_pb05 = &pb05;
    pb05.setX( 10 );
    pb05.setY( pb04.y() + 10 );
    pb05.setWidth( 120 );
    pb05.setMargin( 5 );
    pb05.setTextMargin( 5 );
    pb05.setCaption("Buffer Line");
    //pb04.deactivate();
    pb05.show();

    PB06 pb06( &control_screen );
    p_pb06 = &pb06;
    pb06.setX( pb05.x() + pb05.width() );
    pb06.setY( pb05.y() );
    pb06.setWidth( 80 );
    pb06.setMargin( 5 );
    pb06.setTextMargin( 5 );
    pb06.setCaption("Clear");
    //pb04.deactivate();
    pb06.show();

    PB07 pb07( &control_screen );
    p_pb07 = &pb07;
    pb07.setX( pb06.x() + pb06.width() );
    pb07.setY( pb06.y() );
    pb07.setWidth( 160 );
    pb07.setMargin( 5 );
    pb07.setTextMargin( 5 );
    pb07.setCaption("Save CellTree.dat");
    //pb04.deactivate();
    pb07.show();

    RB00 rb00( &control_screen );
    rb00.setX( pb05.x() );
    rb00.setY( pb05.y() - pb05.height() );
    rb00.setWidth( 100 );
    rb00.setMargin( 15 );
    rb00.setCaption( "Principal" );
    rb00.setState( "true" );

    RB01 rb01( &control_screen );
    rb01.setX( rb00.x() + rb00.width() );
    rb01.setY( rb00.y() );
    rb01.setWidth( 130 );
    rb01.setMargin( 15 );
    rb01.setCaption( "Intermediate" );

    RB02 rb02( &control_screen );
    rb02.setX( rb01.x() + rb01.width() );
    rb02.setY( rb01.y() );
    rb02.setWidth( 100 );
    rb02.setMargin( 15 );
    rb02.setCaption( "Minor" );

    kvs::glut::RadioButtonGroup rbg0( &control_screen );
    rbg0.add( &rb00 );
    rbg0.add( &rb01 );
    rbg0.add( &rb02 );
    rbg0.show();

    CB00 cb00( &control_screen );
    cb00.setX( rb00.x() + rb00.width() );
    cb00.setY( rb00.y() - rb00.height() );
    cb00.setWidth( 100 );
    cb00.setMargin( 25 );
    cb00.setCaption( "Cache" );
    cb00.setState( true );
    cb00.show();

    kvs::ColorMap cmap;
    cmap.setRange( -127, 127 );
    cmap.create();

    kvs::glut::LegendBar legend_bar( &main_screen );
    legend_bar.setColorMap( cmap );
    legend_bar.show();

    return( app.run() );

}
