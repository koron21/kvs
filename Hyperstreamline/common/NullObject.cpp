//
//  NullObject.cpp
//  
//
//  Created by Jun Nishimura on 11/4/10.
//  Copyright 2011 Jun Nishimura. All rights reserved.
//

#include "NullObject.h"

namespace kvs
{

NullObject::NullObject( void ) :
    kvs::ObjectBase()
{
}

NullObject::NullObject( const kvs::ObjectBase* object ) :
    kvs::ObjectBase()
{
    const kvs::Vector3f min_coord = object->minObjectCoord();
    const kvs::Vector3f max_coord = object->maxObjectCoord();

    BaseClass::setMinMaxObjectCoords( min_coord, max_coord );
    BaseClass::setMinMaxExternalCoords( min_coord, max_coord );
}

NullObject::NullObject( const kvs::Vector3f& min_coord, const kvs::Vector3f& max_coord ) :
    kvs::ObjectBase()
{
    BaseClass::setMinMaxObjectCoords( min_coord, max_coord );
    BaseClass::setMinMaxExternalCoords( min_coord, max_coord );
}

NullObject::~NullObject( void )
{
}

kvs::NullObject* NullObject::DownCast( kvs::ObjectBase* object )
{
    const kvs::ObjectBase::ObjectType type = object->objectType();
    if ( type != kvs::ObjectBase::Geometry )
    {
        kvsMessageError( "Input object is not a geometry object." );
        return( NULL );
    }

    kvs::NullObject* null = static_cast<kvs::NullObject*>( object );

    return( null );
}

const kvs::NullObject* NullObject::DownCast( const kvs::ObjectBase* object )
{
    return( NullObject::DownCast( const_cast<kvs::ObjectBase*>( object ) ) );
}

void NullObject::shallowCopy( const NullObject& object )
{
    BaseClass::operator=( object );
}

void NullObject::deepCopy( const NullObject& object )
{
    BaseClass::operator=( object );
}

const kvs::ObjectBase::ObjectType NullObject::objectType( void ) const
{
    return( kvs::ObjectBase::Geometry );
}

}
