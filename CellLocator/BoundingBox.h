#ifndef KVS__BOUNDING_BOX_H_INCLUDE
#define KVS__BOUNDING_BOX_H_INCLUDE

#include <kvs/UnstructuredVolumeObject>

namespace kvs
{
    class BoundingBox
    {
    public:
        BoundingBox();
        ~BoundingBox();
        BoundingBox(const kvs::UnstructuredVolumeObject* object, size_t cindex);		///<	create a bounding box based on a specific object and cell index 
    public:
        const float* bounds() const;
    private:
        float m_bounds[6];
    };
}

#endif
