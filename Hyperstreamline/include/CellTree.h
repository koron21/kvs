#ifndef KVS__CELL_TREE_H_INCLUDE
#define KVS__CELL_TREE_H_INCLUDE

#include "BoundingBox.h"
#include <kvs/UnstructuredVolumeObject>
#include <kvs/Thread>
#include <kvs/BitArray>

namespace kvs
{
    struct stack
    {
        unsigned int    m_stack[32];
        unsigned int*   m_p_stack;

        stack()
        {
            m_p_stack = m_stack;
        }

        ~stack()
        {
        }

        unsigned int size()
        {
            return ( m_p_stack - m_stack );
        }

        void push( const unsigned int num )
        {
            if ( size() > 31 )
            {
                std::cout << "celltree.h: supportive stack overflow!" << std::endl;
                exit(1);
            }

            *(m_p_stack++) = num;
        }

        const unsigned int pop()
        {
            if ( m_p_stack == m_stack )
            {
                std::cout << "celltree.h: supportvie stack bottom reached!" << std::endl;
                exit(0);
            }

            return *(--m_p_stack);
        }

        const unsigned int top()
        {
            if ( m_p_stack == m_stack )
            {
                std::cout << "celltree.h: supportvie stack bottom reached!" << std::endl;
                exit(0);
            }

            return *m_p_stack;
        }

        bool has( const unsigned int num ) const
        {
            unsigned int* p_stack = m_p_stack;

            while ( p_stack-- != m_stack )
            {
                if (*p_stack == num)
                    return true;
            }

            return false;
        }

    };

    struct bucket // only stores current max/min and time of add()
    {
        float        min;
        float        max;
        unsigned int cnt;    //counter of add()
    
        bucket()
        {
            cnt = 0;
            min =  std::numeric_limits<float>::max();
            max = -std::numeric_limits<float>::max();
        }
    
        void add( const float _min, const float _max )
        {
            ++cnt;
        
            if( _min < min )
                min = _min;
            
            if( _max > max )
                max = _max;
        }
    };

    struct per_cell 
    {
        float        min[3];
        float        max[3];
        unsigned int ind;
    };

    struct center_order
    {
        unsigned int d;
    
        center_order( unsigned int _d ) : 
            d(_d)
        {
        }

        bool operator()( const per_cell& pc0, const per_cell& pc1 )
        {
            return (pc0.min[d] + pc0.max[d]) < (pc1.min[d] + pc1.max[d]);
        }
    };

    struct left_predicate
    {
        unsigned int       d;
        float              p;

        left_predicate()
        {
            d = 0;
            p = 0;
        }
    
        left_predicate( unsigned int _d, float _p ) : 
            d(_d), p(2.0f*_p)
        {
        }
   
        bool operator()( const per_cell& pc )
        {
            return (pc.min[d] + pc.max[d]) < p;
        }
    };

    void find_min_max( const per_cell* begin, const per_cell* end, float* min, float* max );
    
    void find_min_d( const per_cell* begin, const per_cell* end, unsigned int d, float& min );

    void find_max_d( const per_cell* begin, const per_cell* end, unsigned int d, float& max );

class CellTree
{
public:
    CellTree(){}
    ~CellTree()
    {
        nodes.clear();
        leaves.clear();
    }

public:
    struct node
    {
        unsigned int index;

        union {
            struct {
                float lm;
                float rm;
            };

            struct {
                unsigned int sz;
                unsigned int st;
            };
        };

        // index = dim | m_nodes.size() * 4
        // so index minimun is 000 | 100 = 101 = 4
        void make_node( unsigned int left, unsigned int d, float b[2] )
        {
            index = (d & 3) | (left << 2);
            lm = b[0];
            rm = b[1];
        }

        void set_children( unsigned int left )
        {
            index = dim() | (left << 2);
        }

        bool is_node() const
        {
            return (index & 3) != 3;
        }

        unsigned int left() const
        {
            return (index >> 2);
        }

        unsigned int right() const
        {
            return (index >> 2) + 1;
        }

        unsigned int dim() const
        {
            return index & 3;
        }

        const float& lmax() const
        {
            return lm;
        }

        const float& rmin() const
        {
            return rm;
        }


        void make_leaf( unsigned int start, unsigned int size )
        {
            index = 3;
            sz = size;
            st = start;
        }

        bool is_leaf() const
        {
            return index == 3;
        }

        unsigned int start() const
        {
            return st;
        }

        unsigned int size() const
        {
            return sz;
        }
    };
      
    std::vector<node>         nodes;
    std::vector<unsigned int> leaves;

    struct pre_traversal
    {
        const CellTree&        m_ct;
        unsigned int           m_stack[32];
        unsigned int*          m_sp;
        const float*           m_pos;

        pre_traversal( const CellTree& ct, const float* pos ) :
            m_ct(ct), m_pos(pos)
        {
            m_stack[0] = 0;             // initialize stack pop to zero
            m_sp = m_stack + 1;         // initialize stack pointer 
        }

        const CellTree::node* next()
        {
            while( true )
            {
                 if( m_sp == m_stack )   //if stack is empty
                    return 0;

                // &m_ct.nodes.front() is the address of the root node
                // move the stack pointer backward by 1, and take it's dereference
                const CellTree::node* n = &m_ct.nodes.front() + *(--m_sp);
                
                if( n->is_leaf() )  
                    return n;

                const float p = m_pos[n->dim()];        // the value corresponding to split dim of the node
                const unsigned int left = n->left();    // get the VECTOR STORAGE index of the left child

                bool l = p <= n->lmax();
                bool r = p > n->rmin();

                if( l && r )
                {
                    if( n->lmax()-p < p-n->rmin() )
                    {
                        *(m_sp++) = left;
                        *(m_sp++) = left+1;
                    }
                    else
                    {
                        *(m_sp++) = left+1;
                        *(m_sp++) = left;
                    }
                }
                else if( l )
                    *(m_sp++) = left;
                else if( r )
                    *(m_sp++) = left+1;
            }
        }
    };

    struct pre_traversal_cached
    {
        const CellTree&         m_ct;
        unsigned int            m_stack[32];
        unsigned int*           m_sp;
        const float*            m_pos;

        // dangerous, be sure hint_stack and hint_sp are well initialized
        pre_traversal_cached( const CellTree& ct, const float* pos, const unsigned int previous_index ) :
            m_ct(ct), m_pos(pos)
        {
            m_stack[0] = 0;                 // initialize stack pointer 
            if ( previous_index != 0 )
            {  
                m_stack[1] = previous_index;
                m_sp = m_stack + 2;
            }
            else
            {
                m_sp = m_stack + 1;
            }
        }

        const CellTree::node* next()
        {
            while( true )
            {
                if( m_sp == m_stack )   //if stack is empty
                    return 0;

                // &m_ct.nodes.front() is the address of the root node
                // move the stack pointer backward by 1, and take it's dereference
                const CellTree::node* n = &m_ct.nodes.front() + *(--m_sp);
                
                if( n->is_leaf() )  
                    return n;
                
                const float p = m_pos[n->dim()];        // the value corresponding to split dim of the node
                const unsigned int left = n->left();    // get the VECTOR STORAGE index of the left child

                bool l = p <= n->lmax();
                bool r = p > n->rmin();

                if( l && r )
                {
                    if( n->lmax()-p < p-n->rmin() )
                    {
                        *(m_sp++) = left;
                        *(m_sp++) = left+1;
                    }
                    else
                    {
                        *(m_sp++) = left+1;
                        *(m_sp++) = left;
                    }
                }
                else if( l )
                    *(m_sp++) = left;
                else if( r )
                    *(m_sp++) = left+1;
            }


        }

        const unsigned int* stack()
        {
            return m_stack;
        }

        const unsigned int* sp()
        {
            return m_sp;
        }

    };

    struct in_traversal_cached
    {
        const CellTree&         m_ct;
        unsigned int            m_stack[64];
        unsigned int*           m_sp;
        const float*            m_pos;
        stack                   m_lrstack;
  
        in_traversal_cached( const CellTree& ct, const float* pos, unsigned int hint_stack[64], 
            unsigned int* hint_sp ):
            m_ct( ct ), m_pos( pos )
        {
            memcpy( m_stack, hint_stack, 256 );

            int n = hint_sp - hint_stack;     // initialize stack pointer 
            m_sp = m_stack + n;
         }

        const CellTree::node* next()
        {
            while( true )
            {

                if( m_sp == m_stack )   //if returned to the root
                    return 0;

                // &m_ct.nodes.front() is the address of the root node
                // move the stack pointer backward by 1, and take it's dereference
                const CellTree::node* n = &m_ct.nodes.front() + *(--m_sp);
                
                if( n->is_leaf() ) 
                {
                    return n;
                }
                
                const float p = m_pos[n->dim()];        // the value corresponding to split dim of the node
                const unsigned int left = n->left();    // get the VECTOR STORAGE index of the left child

                bool l = p <= n->lmax();
                bool r = p > n->rmin();

                if( l && r )
                {
                    if ( !m_lrstack.has( *m_sp ) )
                    {
                        m_lrstack.push( *m_sp );
                    }
                    else // if already registered in lr_stack
                    {
                        *( m_sp + 1 ) = -1;
                        *( m_sp + 2 ) = -1;
                        continue;
                    }

                    if( n->lmax()-p < p-n->rmin() )
                    {
                        m_sp++;
                        *(m_sp++) = left;
                        *(m_sp++) = left + 1;
                    }
                    else
                    {
                        m_sp++;
                        *(m_sp++) = left + 1;
                        *(m_sp++) = left;
                    }
                }
                else if( l )
                {
                    if ( *( m_sp + 1 ) == left )
                    {
                        *( m_sp + 1 ) = -1;
                        *( m_sp + 2 ) = -1;
                        continue;
                    }
                    m_sp++;
                    *(m_sp++) = left;
                }
                else if( r )
                {
                    if ( *( m_sp + 1 ) == left+1 )
                    {
                        *( m_sp + 1 ) = -1;
                        *( m_sp + 2 ) = -1;
                        continue;      
                    }
                    m_sp++;
                    *(m_sp++) = left+1;
                }
            }
 
        }
    };

    unsigned int height( node& node )
    {
        if ( node.index == 3 )
            return 0;
        else
            return ( 1 + std::max( height( nodes[ node.left() ] ), height( nodes[ node.right() ] ) ) );

    }
 
};

class SplitThread : public kvs::Thread
{
public:
    SplitThread(){}
    ~SplitThread(){}

public:
    void init(
        unsigned int leafsize,
        std::vector<CellTree::node>* p_nodes,
        per_cell* pc,
        unsigned int index,
        float min[3],
        float max[3] );

    const bool check();
    void run();
    void split( unsigned int index, float min[3], float max[3] );

private:
    unsigned int                    m_leafsize;
    std::vector<CellTree::node>*    m_nodes;
    unsigned int                    m_index;
    float                           m_min[3];
    float                           m_max[3];
    per_cell*                       m_pc;
};

class CellTreeBuilder
{

public:  
    CellTreeBuilder();
    ~CellTreeBuilder();

    void setParallel();
    void build( CellTree& ct, const kvs::UnstructuredVolumeObject* ds );
    void split( unsigned int index, float min[3], float max[3] );

private:
    bool                                m_parallel;
    unsigned int                        m_leafsize;
    std::vector<kvs::CellTree::node>    m_nodes;
    std::vector<kvs::CellTree::node>    m_nodes1;
    std::vector<kvs::CellTree::node>    m_nodes2;
    kvs::SplitThread                    m_thread[2];
    per_cell*                           m_pc;
    per_cell*                           m_pc1;
    per_cell*                           m_pc2;

};
}


#endif
