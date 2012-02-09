/*
 *  ControlScreen.h
 *  
 *
 *  Created by njun on 11/10/11.
 *  Copyright 2011 Jun Nishimura. All rights reserved.
 *
 */

#ifndef KVS__CONTROL_SCREEN_H_INCLUDE
#define KVS__CONTROL_SCREEN_H_INCLUDE

#include <kvs/glut/Application>
#include <kvs/glut/Screen>

#include "CubicPointObject.h"
#include "StochasticRenderer.h"
#include "StochasticPointRenderer.h"

namespace kvs
{

class ControlScreen : public kvs::glut::Screen
{

    // Class name.
    kvsClassName( kvs::ControlScreen );

    // Module information.
    kvsModuleBaseClass( kvs::glut::Screen );

protected:

    // Reference only.
    kvs::glut::Screen*                  m_screen;
    kvs::CubicPointObject*              m_point;
    kvs::glew::StochasticRenderer*      m_renderer;
    kvs::glew::StochasticPointRenderer* m_point_renderer;

public:

    ControlScreen( kvs::glut::Application* app );

    virtual ~ControlScreen( void );

public:

    void initialize( void );

    void attachMainScreen( 
            kvs::glut::Screen* screen, 
            kvs::CubicPointObject* point, 
            kvs::glew::StochasticRenderer* renderer,
            kvs::glew::StochasticPointRenderer* point_renderer );

public:

    void mouseMoveEvent( kvs::MouseEvent* event );

};

}

#endif
