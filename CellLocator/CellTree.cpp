#include "CellTree.h"

#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>


namespace
{
// of all those cells which pointed by begin and end
// find the overall min max of each dimension, ie, 
// to find the bounding box defined by those cells
// and record min & max
void find_min_max( const per_cell* begin, const per_cell* end, float* min, float* max )
{
    if( begin == end )
        return;
        
    for( unsigned int d=0; d<3; ++d )
    {
        min[d] = begin->min[d];
        max[d] = begin->max[d];
    }
    
    while( ++begin != end )
    {
        for( unsigned int d=0; d<3; ++d )
        {
            if( begin->min[d] < min[d] )    min[d] = begin->min[d];
            if( begin->max[d] > max[d] )    max[d] = begin->max[d];
        }
    }
}

// of all those cells between begin and end
// find the minimum value over a certain dimension
void find_min_d( const per_cell* begin, const per_cell* end, unsigned int d, float& min )
{
    min = begin->min[d];
    
    while( ++begin != end )
        if( begin->min[d] < min )    
            min = begin->min[d];
}

// of all those cells between begin and end
// find the minimum value over a certain dimension
void find_max_d( const per_cell* begin, const per_cell* end, unsigned int d, float& max )
{
    max = begin->max[d];
    
    while( ++begin != end )
        if( begin->max[d] > max )    
            max = begin->max[d];
}
}


namespace kvs
{

CellTreeBuilder::CellTreeBuilder()
{
	m_parallel = false;
    m_leafsize = 8;
}

CellTreeBuilder::~CellTreeBuilder()
{
	m_nodes.clear();
	m_nodes1.clear();
	m_nodes2.clear();
	delete[] m_pc;
}

void CellTreeBuilder::setParallel()
{
	m_parallel = true;
}

void CellTreeBuilder::build( CellTree& ct, const kvs::UnstructuredVolumeObject* ds )
{
    const size_t ncells  = ds->ncells();
    assert( ncells <= std::numeric_limits<unsigned int>::max() );

    m_pc = new per_cell[ncells];

    // max[3], min[3] are the bounding box of the whole data
    float min[3] = { 
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max()
    };

    float max[3] = { 
        -std::numeric_limits<float>::max(),
        -std::numeric_limits<float>::max(),
        -std::numeric_limits<float>::max(),
    };
    
    // m_pc coressponds to each cell
    for( unsigned int i=0; i<ncells; ++i )
    {
        m_pc[i].ind = i;

        //double bounds[6];
        //ds.GetCellBounds( i, bounds );
        //bounds[6] is like: xmin, xmax, ymin, ymax, zmin, zmax
        kvs::BoundingBox bd( ds, i );
        
        const float* bounds = bd.bounds();
        
        /*debug for bounding box
        for( int i=0; i<6; i++ )
            std::cout << " " << bounds[i] ;
        std::cout << std::endl;
        */

        for( int d=0; d<3; ++d )
        {
            m_pc[i].min[d] = bounds[2*d+0];
            m_pc[i].max[d] = bounds[2*d+1];

            if( m_pc[i].min[d] < min[d] )
                min[d] = m_pc[i].min[d];

            if( m_pc[i].min[d] > max[d] )
                max[d] = m_pc[i].max[d];
        }
    }
            
    /* debug for the overall bounding box */
    std::cout << "max external coords found by kvs::BoundingBox:\t" << max[0] << " " << max[1] << " "  << max[2] << " "  << std::endl;
    std::cout << "min external coords found by kvs::BoundingBox:\t" << min[0] << " " << min[1] << " "  << min[2] << " "  << std::endl;


    if ( !m_parallel )
    {
        CellTree::node root;
        root.make_leaf( 0, ncells );	//set index = 3 start = 0 size = ncells
        m_nodes.push_back( root );		

        split( 0, min, max );			//do the split for all the nodes
    
        // rearrange the node's order and copy it to cell tree
        // since the split algorithm do in a pre-order way, a tree such as
        //            0
        //       |----|-----|
        //      1           2
        //   |--|--|     |--|--|
        //   3     4     5     6
        //       |-|-|       |-|-|
        //       7   8       9   10
        // will be stored as an arry like
        // 0 1 2 3 4 7 8 5 6 9 10
        // rearrange the elements such that the ones copied into celltree is like
        // 0 1 2 3 4 5 6 7 8 9 10
        // i.e, in a beautiful up-to-down order
        // and of course, reset their (the node's) indices so they always point to the correct left child
        ct.nodes.resize( m_nodes.size() );
		std::cout << "nnodes: " << m_nodes.size() << std::endl;
        ct.nodes[0] = m_nodes[0];

        std::vector<CellTree::node>::iterator ni = ct.nodes.begin();
        std::vector<CellTree::node>::iterator nn = ct.nodes.begin()+1;

        for( ; ni!=ct.nodes.end(); ++ni )
        {
            if( ni->is_leaf() )
                continue;
        
            *(nn++) = m_nodes[ni->left()];
            *(nn++) = m_nodes[ni->right()];

            ni->set_children( nn-ct.nodes.begin()-2 );
        }

        // --- final 
    
        ct.leaves.resize( ncells );

        for( unsigned int i=0; i<ncells; ++i )
            ct.leaves[i] = m_pc[i].ind;  
    }
    else
    {
        CellTree::node root;
        root.make_leaf( 0, ncells );	//set index = 3 start = 0 size = ncells

        // seperate m_pc into 2 branches
        unsigned int start = 0;
        unsigned int size  = ncells;
    
        per_cell* begin = m_pc + start;
        per_cell* end   = m_pc + start + size;
        per_cell* mid   = begin;

        const unsigned int nbuckets = 6;
        // bounding box calculated by the max, min passed to the function
        const float ext[3] = { max[0]-min[0], max[1]-min[1], max[2]-min[2] }; 
        // iex[3] = { 6/lengthx, 6/lengthy, 6/lengthz }
        const float iext[3] = { nbuckets/ext[0], nbuckets/ext[1], nbuckets/ext[2] };

        bucket b[3][nbuckets];     //3 dimensions, 6 buckets each

        // This part loops though all the cells in range [ begin, end )
        // in x, y, z dimension for only once
        // then counts the number of cells in each buckets, and also,
        // the max and min values of them
        for( const per_cell* pc=begin; pc!=end; ++pc )
        {
            for( unsigned int d=0; d<3; ++d )
            {
                float cen = (pc->min[d] + pc->max[d])/2.0f;	//center of a certain dimension
                // (distance from center to min) / (distance from max to min) * nbuckets, of a certain dimension
                int   ind = (int)( (cen-min[d])*iext[d] );  

                if( ind<0 )	// how can ind < 0 ??? 
                    ind = 0;

                if( ind>=nbuckets ) // how can ind > nbuckets??
                    ind = nbuckets-1;

                b[d][ind].add( pc->min[d], pc->max[d] );
            }
        }
    
        float cost = std::numeric_limits<float>::max();
        float plane;
        unsigned int dim;


        // This part loops through every buckets(6) in all the dimension(3)
        // to determine the best spliting plane
        // bucket1   bucket2   bucket3   bucket4   bucket5   bucket6
        //         |        |         |         |         |
        // each plane is evaluated for cost which is determined by minimizing 
        // left volume * left nnodes + right volume * right nnodes

        // for the fist split, a proper dimension is chosen for the least cost
        for( unsigned int d=0; d<3; ++d )    
        {
            unsigned int sum = 0;
        
            for( unsigned int n=0; n<nbuckets-1; ++n )
            {
                float lmax = -std::numeric_limits<float>::max();
                float rmin =  std::numeric_limits<float>::max();

                for( unsigned int m=0; m<=n; ++m )
                    if( b[d][m].max > lmax )
                        lmax = b[d][m].max;
            
                for( unsigned int m=n+1; m<nbuckets; ++m )
                    if( b[d][m].min < rmin )
                        rmin = b[d][m].min;
            
                sum += b[d][n].cnt;
            
                float lvol = (lmax-min[d])/ext[d];			//left volume
                float rvol = (max[d]-rmin)/ext[d];			//right volume
            
                float c = lvol*sum + rvol*(size-sum);
            
                if( sum > 0 && sum < size && c < cost )
                {
                    cost    = c;
                    dim     = d;
                    plane   = min[d] + (n+1)/iext[d];
                }
            }
        }


        // by using STL algorithm partition,
        // all the cell( containing bounding box, cell id)  is divided to left and right,
        // according to the seperating plane and dimension

        if( cost != std::numeric_limits<float>::max() )		// how can that be ...
            mid = std::partition( begin, end, left_predicate( dim, plane ) );

        // fallback
        // in the extreme case, use split-median to determin the mid
        if( mid == begin || mid == end )
        {
            // the lagest element between ext ext+1 and ext+2
            // note: ext is float type pointer -> 4bytes
            //       dim is uint16 type        -> 4bytes
            dim = std::max_element( ext, ext+3 ) - ext;

            mid = begin + (end-begin)/2;
            std::nth_element( begin, mid, end, center_order( dim ) );
        }

        float lmin[3], lmax[3], rmin[3], rmax[3];

        // find each part's left right bounds
        find_min_max( begin, mid, lmin, lmax );
        find_min_max( mid,   end, rmin, rmax );

        // and choose the one along the splitting dimension
        float clip[2] = { lmax[dim], rmin[dim] };
        
        // assign each to either thread
        unsigned int size1 = mid - m_pc;
        unsigned int start1 = 0;
        unsigned int size2 = end - mid;
        unsigned int start2 = 0;
         
        CellTree::node root1;
        root1.make_leaf( start1, size1 );	//set index = 3 start = 0 size = ncells
        m_nodes1.push_back( root1 );

        CellTree::node root2;
        root2.make_leaf( start2, size2 );	//set index = 3 start = 0 size = ncells
        m_nodes2.push_back( root2 );

        m_pc1 = m_pc;
        m_pc2 = mid;

        m_thread[0].init( m_leafsize, &m_nodes1, m_pc1, 0, lmin, lmax );
	    m_thread[1].init( m_leafsize, &m_nodes2, m_pc2, 0, rmin, rmax );
	
	    m_thread[0].start();
	    m_thread[1].start();

	    m_thread[0].wait();
	    m_thread[1].wait();

        // merge data into celltree
		// size = size_of_tree1 + size_of_tree2 + root
        const unsigned int node_size = m_nodes1.size() + m_nodes2.size() + 1;
        ct.nodes.resize( node_size );
		std::cout << "nnodes: " << node_size << std::endl;

        ct.nodes[0] = root;
        ct.nodes[0].make_node( 1, dim, clip );

        ct.nodes[1] = m_nodes1[0];
        ct.nodes[2] = m_nodes2[0];

        kvs::BitArray mask;
        mask.allocate( node_size );
        mask.reset(1);
        mask.set(2);

        // ni points to root's left child initially
        std::vector<CellTree::node>::iterator ni = ct.nodes.begin()+1;  
        // nn points to toot's left child's left child initially
        std::vector<CellTree::node>::iterator nn = ct.nodes.begin()+3;

        for( unsigned int i = 1; ni!=ct.nodes.end(); ++ni, i++ )
        {
            if( ni->is_leaf() )
                continue;
            if( !mask[i] )
            {
                mask.reset( nn-ct.nodes.begin() );
                *(nn++) = m_nodes1[ni->left()];
                
                mask.reset( nn-ct.nodes.begin() );
                *(nn++) = m_nodes1[ni->right()];
            }
            else
            {
                mask.set( nn-ct.nodes.begin() );
				unsigned int left = ni->left();
				if ( m_nodes2[left].is_leaf() )
					m_nodes2[left].st += size1;
                *(nn++) = m_nodes2[left];

                mask.set( nn-ct.nodes.begin() );
				if ( m_nodes2[left+1].is_leaf() )
					m_nodes2[left+1].st += size1;
                *(nn++) = m_nodes2[left+1];
            }
            ni->set_children( nn-ct.nodes.begin()-2 );
        }

        // --- final 
    
        ct.leaves.resize( ncells );

        for( unsigned int i=0; i<ncells; ++i )
            ct.leaves[i] = m_pc[i].ind;   
    }
#ifdef DEBUG
    std::cout << "nodes array: \n";
    for ( unsigned int i = 0; i < 20; i ++ )
    {
        std::cout << "dim    " << ( ct.nodes[i].index & 3 ) << "\t" 
                  << "left   " << ( ct.nodes[i].index >> 2 ) << "\t"
                  << "start  " << ct.nodes[i].st  << "\t"
                  << "size   " << ct.nodes[i].sz  << "\t"
                  << "lmax   " << ct.nodes[i].lm  << "\t"
                  << "rmin   " << ct.nodes[i].rm << std::endl;
    }
    std::cout << "leaves array: \n";
    for ( unsigned int i = 0; i < 20; i ++ )
    {
        std::cout << "index: " << ct.leaves[i] << std::endl;
    }
#endif


}

void CellTreeBuilder::split( unsigned int index, float min[3], float max[3] )
{
    unsigned int start = m_nodes[index].start();
    unsigned int size  = m_nodes[index].size();
    
    if( size < m_leafsize )		//if size is less than the maxium bucket size, don't do spliting any more
        return;

    per_cell* begin = m_pc + start;
    per_cell* end   = m_pc + start + size;
    per_cell* mid   = begin;
    
    const unsigned int nbuckets  = 6;

    // bounding box calculated by the max, min passed to the function
    const float ext[3] = { max[0]-min[0], max[1]-min[1], max[2]-min[2] }; 
    // iex[3] = { 6/lengthx, 6/lengthy, 6/lengthz }
    const float iext[3] = { nbuckets/ext[0], nbuckets/ext[1], nbuckets/ext[2] };

    bucket b[3][nbuckets]; //3 dimensions, 6 buckets each

    // This part loops though all the cells in range [ begin, end )
    // in x, y, z dimension for only once
    // then counts the number of cells in each buckets, and also,
    // the max and min values of them
    for( const per_cell* pc=begin; pc!=end; ++pc )
    {
        for( unsigned int d=0; d<3; ++d )
        {
            float cen = (pc->min[d] + pc->max[d])/2.0f;	//center of a certain dimension
            // (distance from center to min) / (distance from max to min) * nbuckets, of a certain dimension
            int   ind = (int)( (cen-min[d])*iext[d] );  

            if( ind<0 )	// how can ind < 0 ??? 
                ind = 0;

            if( ind>=nbuckets ) // how can ind > nbuckets??
                ind = nbuckets-1;

            b[d][ind].add( pc->min[d], pc->max[d] );
        }
    }
    
    float cost = std::numeric_limits<float>::max();
    float plane;
    unsigned int dim;


    // This part loops through every buckets(6) in all the dimension(3)
    // to determine the best spliting plane
    // bucket1   bucket2   bucket3   bucket4   bucket5   bucket6
    //         |        |         |         |         |
    // each plane is evaluated for cost which is determined by minimizing 
    // left volume * left nnodes + right volume * right nnodes

    // for the fist split, a proper dimension is chosen for the least cost
    for( unsigned int d=0; d<3; ++d )    
    {
        unsigned int sum = 0;
        
        for( unsigned int n=0; n<nbuckets-1; ++n )
        {
            float lmax = -std::numeric_limits<float>::max();
            float rmin =  std::numeric_limits<float>::max();

            for( unsigned int m=0; m<=n; ++m )
                if( b[d][m].max > lmax )
                    lmax = b[d][m].max;
            
            for( unsigned int m=n+1; m<nbuckets; ++m )
                if( b[d][m].min < rmin )
                    rmin = b[d][m].min;
            
            sum += b[d][n].cnt;
            
            float lvol = (lmax-min[d])/ext[d];			//left volume
            float rvol = (max[d]-rmin)/ext[d];			//right volume
            
            float c = lvol*sum + rvol*(size-sum);
            
            if( sum > 0 && sum < size && c < cost )
            {
                cost    = c;
                dim     = d;
                plane   = min[d] + (n+1)/iext[d];
            }
        }
    }


    // by using STL algorithm partition,
    // all the cell( containing bounding box, cell id)  is divided to left and right,
    // according to the seperating plane and dimension

    if( cost != std::numeric_limits<float>::max() )		// how can that be ...
        mid = std::partition( begin, end, left_predicate( dim, plane ) );

    // fallback
    // in the extreme case, use split-median to determin the mid
    if( mid == begin || mid == end )
    {
        // the lagest element between ext ext+1 and ext+2
        // note: ext is float type pointer -> 4bytes
        //       dim is uint16 type        -> 4bytes
        dim = std::max_element( ext, ext+3 ) - ext;

        mid = begin + (end-begin)/2;
        std::nth_element( begin, mid, end, center_order( dim ) );
    }

    float lmin[3], lmax[3], rmin[3], rmax[3];

    // find each part's left right bounds
    find_min_max( begin, mid, lmin, lmax );
    find_min_max( mid,   end, rmin, rmax );

    // and choose the one along the splitting dimension
    float clip[2] = { lmax[dim], rmin[dim] };

    // making two child leaf nodes, which is consecutive in memory address
    kvs::CellTree::node child[2];
    child[0].make_leaf( begin - m_pc, mid-begin );
    child[1].make_leaf( mid   - m_pc, end-mid );
    
    // then, make the originally leaf node m_nodes[index] an ordinary node
    // index is node.size() + 00 | dim
    // insert the two child nodes at the end of current container 
    m_nodes[index].make_node( m_nodes.size(), dim, clip );
    m_nodes.insert( m_nodes.end(), child, child+2 );

	// traverse to the left brunch
	split( m_nodes[index].left(), lmin, lmax );
	// traverse to the right brunch
	split( m_nodes[index].right(), rmin, rmax );

	
}

void SplitThread::init( 
        unsigned int leafsize,
		std::vector<CellTree::node>* p_nodes,
        per_cell* pc,
		unsigned int index,
		float min[3],
		float max[3] )
{
    m_leafsize = leafsize;
	m_index    = index;
	for ( unsigned int i = 0; i < 3; i ++ )
	{
		m_min[i] = min[i];
		m_max[i] = max[i];
	}
	m_nodes = p_nodes;
    m_pc = pc;
}

const bool SplitThread::check()
{
	return ( &m_nodes != NULL );
}

void SplitThread::run()
{

	split( m_index, m_min, m_max );
}

void SplitThread::split( unsigned int index, float min[3], float max[3] )
{
    std::vector<CellTree::node>& nodes = *m_nodes;
    unsigned int start = nodes[index].start();
    unsigned int size  = nodes[index].size();
    
    if( size < m_leafsize )		//if size is less than the maxium bucket size, don't do spliting any more
        return;

    per_cell* begin = m_pc + start;
    per_cell* end   = m_pc + start + size;
    per_cell* mid   = begin;

    const unsigned int nbuckets = 6;
    // bounding box calculated by the max, min passed to the function
    const float ext[3] = { max[0]-min[0], max[1]-min[1], max[2]-min[2] }; 
    // iex[3] = { 6/lengthx, 6/lengthy, 6/lengthz }
    const float iext[3] = { nbuckets/ext[0], nbuckets/ext[1], nbuckets/ext[2] };

    bucket b[3][nbuckets]; //3 dimensions, 6 buckets each

    // This part loops though all the cells in range [ begin, end )
    // in x, y, z dimension for only once
    // then counts the number of cells in each buckets, and also,
    // the max and min values of them
    for( const per_cell* pc=begin; pc!=end; ++pc )
    {
        for( unsigned int d=0; d<3; ++d )
        {
            float cen = (pc->min[d] + pc->max[d])/2.0f;	//center of a certain dimension
            // (distance from center to min) / (distance from max to min) * nbuckets, of a certain dimension
            int   ind = (int)( (cen-min[d])*iext[d] );  

            if( ind<0 )	// how can ind < 0 ??? 
                ind = 0;

            if( ind>=nbuckets ) // how can ind > nbuckets??
                ind = nbuckets-1;

            b[d][ind].add( pc->min[d], pc->max[d] );
        }
    }
    
    float cost = std::numeric_limits<float>::max();
    float plane;
    unsigned int dim;


    // This part loops through every buckets(6) in all the dimension(3)
    // to determine the best spliting plane
    // bucket1   bucket2   bucket3   bucket4   bucket5   bucket6
    //         |        |         |         |         |
    // each plane is evaluated for cost which is determined by minimizing 
    // left volume * left nodes + right volume * right nnodes

    // for the fist split, a proper dimension is choosed for the least cost
    for( unsigned int d=0; d<3; ++d )    
    {
        unsigned int sum = 0;
        
        for( unsigned int n=0; n<nbuckets-1; ++n )
        {
            float lmax = -std::numeric_limits<float>::max();
            float rmin =  std::numeric_limits<float>::max();

            for( unsigned int m=0; m<=n; ++m )
                if( b[d][m].max > lmax )
                    lmax = b[d][m].max;
            
            for( unsigned int m=n+1; m<nbuckets; ++m )
                if( b[d][m].min < rmin )
                    rmin = b[d][m].min;
            
            sum += b[d][n].cnt;
            
            float lvol = (lmax-min[d])/ext[d];			//left volume
            float rvol = (max[d]-rmin)/ext[d];			//right volume
            
            float c = lvol*sum + rvol*(size-sum);
            
            if( sum > 0 && sum < size && c < cost )
            {
                cost    = c;
                dim     = d;
                plane   = min[d] + (n+1)/iext[d];
            }
        }
    }


    // by using STL algorithm partition,
    // all the cell( containing bounding box, cell id)  is divided to left and right,
    // according to the seperating plane and dimension

    if( cost != std::numeric_limits<float>::max() )		// how can that be ...
        mid = std::partition( begin, end, left_predicate( dim, plane ) );

    // fallback
    // in the extreme case, use split-median to determin the mid
    if( mid == begin || mid == end )
    {
        // the lagest element between ext ext+1 and ext+2
        // note: ext is float type pointer -> 4bytes
        //       dim is uint16 type        -> 4bytes
        dim = std::max_element( ext, ext+3 ) - ext;

        mid = begin + (end-begin)/2;
        std::nth_element( begin, mid, end, center_order( dim ) );
    }

    float lmin[3], lmax[3], rmin[3], rmax[3];

    // find each part's left right bounds
    find_min_max( begin, mid, lmin, lmax );
    find_min_max( mid,   end, rmin, rmax );

    // and choose the one along the splitting dimension
    float clip[2] = { lmax[dim], rmin[dim] };

    // making two child leaf nodes, which is consecutive in memory address
    kvs::CellTree::node child[2];
    child[0].make_leaf( begin - m_pc, mid-begin );
    child[1].make_leaf( mid   - m_pc, end-mid );
    
    // then, make the originally leaf node nodes[index] an ordinary node
    // index is node.size() + 00 | dim
    // insert the two child nodes at the end of current container 
    nodes[index].make_node( nodes.size(), dim, clip );
    nodes.insert( nodes.end(), child, child+2 );

	// traverse to the left brunch
	split( nodes[index].left(), lmin, lmax );
	// traverse to the right brunch
	split( nodes[index].right(), rmin, rmax );
}
}

