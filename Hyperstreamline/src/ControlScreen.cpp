/*
 *  ControlScreen.cpp
 *  
 *
 *  Created by njun on 11/10/11.
 *  Copyright 2011 Jun Nishimura. All rights reserved.
 *
 */

#include "ControlScreen.h"

#include <kvs/Mouse>
#include <kvs/MouseButton>

namespace kvs
{

ControlScreen::ControlScreen( kvs::glut::Application* app ) :
    kvs::glut::Screen( app )
{
    this->initialize();
}

ControlScreen::~ControlScreen( void )
{
}

void ControlScreen::initialize( void )
{
    BaseClass::setTitle( "kvs::ControlScreen" );
    BaseClass::setGeometry( 512, 0, 600, 560 );

    m_screen = NULL;
    m_point = NULL;
    m_renderer = NULL;
    m_point_renderer = NULL;
}

void ControlScreen::attachMainScreen( 
        kvs::glut::Screen* screen, 
        kvs::CubicPointObject* point,
        kvs::glew::StochasticRenderer* renderer,
        kvs::glew::StochasticPointRenderer* point_renderer)
{
    m_screen = screen;
    m_point  = point;
    m_renderer = renderer;
    m_point_renderer = point_renderer;
}


void ControlScreen::mouseMoveEvent( kvs::MouseEvent* event )
{
    if ( event->x() > 20 && event->x() < 330 && event->y() > 20 && event->y() < 400 )
    {
        const kvs::Xform x = m_screen->objectManager()->xform();

        const float* pcoord = m_point->coords().pointer();
        const unsigned int nvertices = m_point->nvertices();
        kvs::ValueArray<float> coords( nvertices * 3 );

        if ( event->button() == kvs::MouseButton::Right )
        {
            m_mouse->setMode( kvs::Mouse::Translation );
            m_mouse->move( event->x(), event->y() );
            kvs::Vector3f translation = m_mouse->translation();
            const kvs::Vector3f normalize = m_screen->objectManager()->normalize();

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
            m_point->setCoords( coords );
        }

        if ( event->button() == kvs::MouseButton::Left )
        {
            m_mouse->setMode( kvs::Mouse::Rotation );
            m_mouse->move( event->x(), event->y() );
            kvs::Matrix33f rotation = m_mouse->rotation().toMatrix();

            for ( unsigned int i = 0; i < nvertices; i ++ )
            {
                kvs::Vector3f coord( pcoord );
                const kvs::Vector3f new_coord = coord * rotation;
                coords[ 3 * i] = new_coord.x();
                coords[ 3 * i + 1] = new_coord.y();
                coords[ 3 * i + 2] = new_coord.z();
                pcoord += 3;
            }
            m_point->setCoords( coords );
        }
        m_renderer->changeObject( m_point, m_point_renderer, false );
    }

    BaseClass::eventHandler()->notify( event );
    BaseClass::redraw();
}



}
