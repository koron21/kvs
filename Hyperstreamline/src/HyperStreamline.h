#ifndef KVS__HYPER_STREAMLINE_H_INCLUDE
#define KVS__HYPER_STREAMLINE_H_INCLUDE

#include<algorithm>

#include <kvs/ClassName>
#include <kvs/Module>
#include <kvs/MapperBase>
#include <kvs/LineObject>
#include <kvs/PointObject>
#include <kvs/StructuredVolumeObject>
#include <kvs/UnstructuredVolumeObject>
#include <kvs/Thread>

#include "CellLocatorBIH.h"

namespace kvs
{

class HyperStreamlineThread;
    
class HyperStreamline : public kvs::MapperBase, public kvs::LineObject
{
    // Class name.
    kvsClassName( kvs::HyperStreamline );

    // Module information.
    kvsModuleBaseClass( kvs::MapperBase );
    kvsModuleSuperClass( kvs::LineObject );

public:

    HyperStreamline( void );

    HyperStreamline(
        const kvs::VolumeObjectBase* volume,
        const kvs::PointObject* seed_points,
        const kvs::TransferFunction& transfer_function );

    virtual ~HyperStreamline( void );

public:

        enum IntegrationMethod
        {
            Euler = 0,
            RungeKutta2nd = 1,
            RungeKutta4th = 2
        };

        enum IntegrationDirection
        {
            ForwardDirection = 0,
            BackwardDirection = 1,
            BothDirections = 2
        };

        enum LocationMethod
        {
            BIH = 0,
            ADJ = 1
        };


public:

    SuperClass* exec( const kvs::ObjectBase* object );

public:
    const bool check_for_inside_volume( const kvs::Vector3f& seed );

public:

    const std::vector<float>& eigenValues() const;

    void setEigenValues( const std::vector<float>& eigen_values );
    void setSeedPoints( const kvs::PointObject* seed_points );
    void setIntegrationMethod( const HyperStreamline::IntegrationMethod method );
    void setIntegrationDirection( const HyperStreamline::IntegrationDirection direction );
    void setIntegrationInterval( const float interval );
    void setVectorLengthThreshold( const float length );
    void setIntegrationTimesThreshold( const size_t times );
    void setEnableBoundaryCondition( const bool enabled );
    void setEnableVectorLengthCondition( const bool enabled );
    void setEnableIntegrationTimesCondition( const bool enabled );
    void setLocationMethod( const HyperStreamline::LocationMethod method, kvs::CellTree* ct = 0, kvs::UnstructuredVolumeObject* volume = 0 );
    void setLocator( const kvs::CellTree* ct = 0, const kvs::UnstructuredVolumeObject* volume = 0 );
    void setDisableCache();
    void setLocatorToInitialized();
    void setGoWithNthEigenVector( const int nth );

    //only use for update
    void calculate_color( const float min_egvalue, const float max_egvalue );
    void calculate_color();

protected:

    void mapping( const kvs::VolumeObjectBase* volume );
    void extract_lines(
        const kvs::VolumeObjectBase* volume );

    const bool calculate_lines(
        std::vector<kvs::Real32>* vertices,
        std::vector<kvs::UInt8>* colors,
        const size_t index );

protected:

    kvs::PointObject*           m_seed_points;
    std::vector<float>          m_eigenvalues;
    IntegrationMethod           m_integration_method;
    IntegrationDirection        m_integration_direction;
    float                       m_integration_interval;
    float                       m_vector_length_threshold;
    size_t                      m_integration_times_threshold;
    bool                        m_enable_boundary_condition;
    bool                        m_enable_vector_length_condition;
    bool                        m_enable_integration_times_condition;
    bool                        m_degenerate;
    kvs::CellLocatorBIH*        m_locator;
    LocationMethod              m_location_method;
    size_t                      m_nthreads;
    HyperStreamlineThread*      m_p_threads;
    bool                        m_locator_initialized;
    bool                        m_enable_cache;
    int                         m_nth_egvector;

};

class HyperStreamlineThread : public kvs::Thread
{
    // Class name.
    kvsClassName( kvs::HyperStreamlineThread );

    // Module information.
    kvsModuleBaseClass( kvs::Thread );


public:

    HyperStreamlineThread();
    ~HyperStreamlineThread();

public:

    void init();
    const bool check();
    void run();

public:

    const kvs::LineObject& line() const;
    const std::vector<float>& eigenValue() const;

    void setVolume( const kvs::VolumeObjectBase* volume );
    void setTransferFunction( const kvs::TransferFunction& transfer_function );
    void setSeedPoint( const kvs::Vector3f seed );
    void setLocator( const kvs::CellTree* ct, const kvs::UnstructuredVolumeObject* volume );
    void setDirection( const kvs::HyperStreamline::IntegrationDirection Direction);
    void setMethod( const kvs::HyperStreamline::IntegrationMethod Method );
    void setGoWithNthEigenVector( const int nth );
    void setDisableCache();

public:

    const bool calculate_one_side(
        std::vector<kvs::Real32>& coords,
        std::vector<kvs::UInt8>& colors,
        const kvs::Vector3f& seed_point,
        const kvs::Vector3f& seed_vector );

    const bool calculate_next_vertex(
        const kvs::Vector3f& current_vertex,
        const kvs::Vector3f& current_direction,
        kvs::Vector3f* next_vertex );

    const bool integrate_by_euler(
        const kvs::Vector3f& current_vertex,
        const kvs::Vector3f& current_direction,
        kvs::Vector3f* next_vertex );

    const bool integrate_by_runge_kutta_2nd(
        const kvs::Vector3f& current_vertex,
        const kvs::Vector3f& current_direction,
        kvs::Vector3f* next_vertex );

    const bool integrate_by_runge_kutta_4th(
        const kvs::Vector3f& current_vertex,
        const kvs::Vector3f& current_direction,
        kvs::Vector3f* next_vertex );

protected:

    const bool check_for_inside_volume( const kvs::Vector3f& point );
    const bool check_for_vector_length( const kvs::Vector3f& direction );
    const bool check_for_integration_times( const size_t times );

    const bool check_for_acceptance( const std::vector<kvs::Real32>& vertices );
    const bool check_for_termination(
        const kvs::Vector3f& current_vertex,
        const kvs::Vector3f& direction,
        const size_t integration_times,
        const kvs::Vector3f& next_vertex );
    const kvs::Vector3f interpolate_vector( const kvs::Vector3f& vertex, const kvs::Vector3f& direction, float& eigenvalue );
    const kvs::Vector3f calculate_vector( const kvs::Vector3f& vertex );
    const kvs::RGBColor calculate_color();

protected:

    std::vector<float>                              m_eigenvalues;
    float                                           m_mises;
    const kvs::VolumeObjectBase*                    m_volume;
    kvs::TransferFunction                           m_transfer_function;
    kvs::Vector3f                                   m_seed_point;
    kvs::LineObject                                 m_line;
    mutable kvs::CellLocatorBIH*                    m_locator;
    kvs::HyperStreamline::IntegrationDirection      m_direction;
    kvs::HyperStreamline::IntegrationMethod         m_method;

    float                                           m_integration_interval;
    float                                           m_vector_length_threshold;
    size_t                                          m_integration_times_threshold;
    bool                                            m_initialized;
    bool                                            m_enable_boundary_condition;
    bool                                            m_enable_vector_length_condition;
    bool                                            m_enable_integration_times_condition;
    bool                                            m_degenerate;
    int                                             m_nth_egvector;
};


} // end of namespace kvs

#endif // KVS__Hyper_STREAMLINE_BASE_H_INCLUDE
