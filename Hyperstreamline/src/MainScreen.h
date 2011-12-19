/*
 *  MainScreen.h
 *  
 *
 *  Created by njun on 11/10/11.
 *  Copyright 2011 Jun Nishimura. All rights reserved.
 *
 */

#ifndef KVS__MAIN_SCREEN_H_INCLUDE
#define KVS__MAIN_SCREEN_H_INCLUDE

#include <kvs/glut/Application>
#include <kvs/glut/Screen>

#include <kvs/CommandLine>
#include <kvs/ExternalFaces>
#include <kvs/PolygonImporter>
#include <kvs/UnstructuredVolumeImporter>

#include "CellLocatorBIH.h"
#include "HyperStreamline.h"
#include "NullObject.h"
#include "PolygonToPolygon.h"
#include "StochasticLineRenderer.h"
#include "StochasticPointRenderer.h"
#include "StochasticPolygonRenderer.h"
#include "StochasticRenderer.h"
#include "StochasticVolumeRenderer.h"

namespace kvs
{

class MainScreen : public kvs::glut::Screen
{

    // Class name.
    kvsClassName( kvs::MainScreen );

    // Module information.
    kvsModuleBaseClass( kvs::glut::Screen );

protected:

    kvs::UnstructuredVolumeObject*      m_volume1;
    kvs::UnstructuredVolumeObject*      m_volume2;
    kvs::CellLocatorBIH*                m_locator;
    kvs::PointObject*                   m_seed_point;
    kvs::PolygonObject*                 m_polygon;
    kvs::HyperStreamline*               m_streamline;
    kvs::TransferFunction               m_tfunc;
    std::vector<unsigned int>           m_buffer;
    kvs::HyperStreamline*               m_buffered_streamline;

    kvs::glew::StochasticRenderer*          m_renderer;
    kvs::glew::StochasticPointRenderer*     m_point_renderer;
    kvs::glew::StochasticLineRenderer*      m_line_renderer;
    kvs::glew::StochasticLineRenderer*      m_buffered_line_renderer;
    kvs::glew::StochasticPolygonRenderer*   m_polygon_renderer;
    kvs::glew::StochasticVolumeRenderer*    m_volume_renderer;


public:

    MainScreen( kvs::glut::Application* app );
    virtual ~MainScreen( void );

public:
    
    void redraw();
    void keyPressEvent( kvs::KeyEvent* event );
    void initialize( kvs::glut::Application* app );
    void clear( void );

public:

    kvs::UnstructuredVolumeObject* volume1();
    kvs::UnstructuredVolumeObject* volume2();
    kvs::glew::StochasticVolumeRenderer* volumeRenderer();
    kvs::glew::StochasticRenderer* renderer();

    kvs::PointObject* seedPoint( void );
    void updateStreamLine( void );

};

}

#endif
