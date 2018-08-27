#include <map>
#include <fstream>
#include "Mesh.h"
#include "MeshIO.h"
#include "DiscreteExteriorCalculus.h"
#include <random>
#include <functional>

static std::default_random_engine re(123);

using namespace std;

namespace DDG
{
    Mesh :: Mesh( void )
    {}

    Mesh :: Mesh( const Mesh& mesh )
    {
        *this = mesh;
    }

    class  HalfEdgeIterCompare { public: bool operator()( const  HalfEdgeIter& i, const  HalfEdgeIter& j ) const { return &*i < &*j; } };
    class HalfEdgeCIterCompare { public: bool operator()( const HalfEdgeCIter& i, const HalfEdgeCIter& j ) const { return &*i < &*j; } };
    class    VertexIterCompare { public: bool operator()( const    VertexIter& i, const    VertexIter& j ) const { return &*i < &*j; } };
    class   VertexCIterCompare { public: bool operator()( const   VertexCIter& i, const   VertexCIter& j ) const { return &*i < &*j; } };
    class      FaceIterCompare { public: bool operator()( const      FaceIter& i, const      FaceIter& j ) const { return &*i < &*j; } };
    class     FaceCIterCompare { public: bool operator()( const     FaceCIter& i, const     FaceCIter& j ) const { return &*i < &*j; } };
    class      EdgeIterCompare { public: bool operator()( const      EdgeIter& i, const      EdgeIter& j ) const { return &*i < &*j; } };
    class     EdgeCIterCompare { public: bool operator()( const     EdgeCIter& i, const     EdgeCIter& j ) const { return &*i < &*j; } };

    const Mesh& Mesh :: operator=( const Mesh& mesh )
    {
        map< HalfEdgeCIter, HalfEdgeIter, HalfEdgeCIterCompare > halfedgeOldToNew;
        map<   VertexCIter,   VertexIter,   VertexCIterCompare >   vertexOldToNew;
        map<     EdgeCIter,     EdgeIter,     EdgeCIterCompare >     edgeOldToNew;
        map<     FaceCIter,     FaceIter,     FaceCIterCompare >     faceOldToNew;

        // copy geometry from the original mesh and create a
        // map from pointers in the original mesh to
        // those in the new mesh
        halfedges.clear();
        for( HalfEdgeCIter he = mesh.halfedges.begin(); he != mesh.halfedges.end(); he++ )
            halfedgeOldToNew[ he ] = halfedges.insert( halfedges.end(), *he );

        vertices.clear();
        for(   VertexCIter  v =  mesh.vertices.begin();  v !=  mesh.vertices.end();  v++ )
            vertexOldToNew[ v  ] =  vertices.insert(  vertices.end(), *v  );

        edges.clear();
        for(     EdgeCIter  e =     mesh.edges.begin();  e !=     mesh.edges.end();  e++ )
            edgeOldToNew[ e  ] =     edges.insert(     edges.end(), *e  );

        faces.clear();
        for(     FaceCIter  f =     mesh.faces.begin();  f !=     mesh.faces.end();  f++ )
            faceOldToNew[ f  ] =     faces.insert(     faces.end(), *f  );

        // "search and replace" old pointers with new ones
        for( HalfEdgeIter he = halfedges.begin(); he != halfedges.end(); he++ )
        {
            he->next   = halfedgeOldToNew[ he->next   ];
            he->flip   = halfedgeOldToNew[ he->flip   ];
            he->vertex =   vertexOldToNew[ he->vertex ];
            he->edge   =     edgeOldToNew[ he->edge   ];
            he->face   =     faceOldToNew[ he->face   ];
        }

        for( VertexIter v = vertices.begin(); v != vertices.end(); v++ ) v->he = halfedgeOldToNew[ v->he ];
        for(   EdgeIter e =    edges.begin(); e !=    edges.end(); e++ ) e->he = halfedgeOldToNew[ e->he ];
        for(   FaceIter f =    faces.begin(); f !=    faces.end(); f++ ) f->he = halfedgeOldToNew[ f->he ];

        return *this;
    }

    int Mesh::read( const string& filename )
    {
        inputFilename = filename;
        ifstream in( filename.c_str() );

        if( !in.is_open() )
        {
            cerr << "Error reading from mesh file " << filename << endl;
            return 1;
        }

        int rval;
        if( !( rval = MeshIO::read( in, *this )))
        {
            indexElements();
            normalize();
        }
        return rval;
    }

    int Mesh::writeMesh( const string& filename ) const
        // reads a mesh from a Wavefront OBJ file; return value is nonzero
        // only if there was an error
    {
        cout << "writing mesh " << filename << endl;
        ofstream out( filename.c_str() );

        if( !out.is_open() )
        {
            cerr << "Error writing to mesh file " << filename << endl;
            return 1;
        }

        MeshIO::writeMesh( out, *this );

        return 0;
    }

    int Mesh::writeEdgeTube( const std::string& filename, const double radius, const int faces ) const
    {
        cout << "writing tube mesh " << filename << endl;
        ofstream out( filename.c_str() );

        if( !out.is_open() )
        {
            cerr << "Error writing to mesh file " << filename << endl;
            return 1;
        }

        MeshIO::writeEdgeTube( out, *this, radius, faces );

        return 0;
    }
    int Mesh::writeEdgeStrips( const std::string& filename) const
    {

        MeshIO::writeEdgeStrips( filename, *this );

        return 0;             
    }

    int Mesh::writeEdgeCSV( const std::string& filename ) const
    {
        cout << "writing edge csv " << filename << endl;
        ofstream out( filename.c_str() );

        if( !out.is_open() )
        {
            cerr << "Error writing to mesh file " << filename << endl;
            return 1;
        }

        MeshIO::writeEdgeCSV( out, *this );

        return 0;
    }

    bool Mesh::reload( void )
    {
        return read( inputFilename );
    }

    void Mesh::normalize( void )
    {
        // compute center of mass
        Vector c( 0., 0., 0. );
        for( VertexCIter v = vertices.begin(); v != vertices.end(); v++ )
        {
            c += v->position;
        }
        c /= (double) vertices.size();

        // translate to origin
        for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
        {
            v->position -= c;
        }

        // rescale such that the mesh sits inside the unit ball
        double rMax = 0.;
        for( VertexCIter v = vertices.begin(); v != vertices.end(); v++ )
        {
            rMax = max( rMax, v->position.norm() );
        }
        for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
        {
            v->position /= rMax;
        }

        updateNormals();
    }

    void Mesh::noise(const double sigma_raw)
    {
        cout << "noise " << sigma_raw << endl;
        const double sigma = meanEdgeLength()*sigma_raw;
        auto dist = std::bind(std::normal_distribution<double>(0, sigma), re);

        for( VertexIter vertex = vertices.begin(); vertex != vertices.end(); vertex++ )
        {
            const Vector random_stuff(dist(), dist(), dist());
            vertex->position += random_stuff;
        }

        updateNormals();
    }

    void Mesh::updateNormals()
    {
        int count = 0;
        for (Face& face : faces)
        {
            face.normal = face.from_vertices_normal();
            if (face.normal.norm()) continue;
            count++;
            face.normal = face.geometric_normal();
        }
        cout << count << "/" << faces.size() << " geometric normals" << endl;
    }

    void Mesh::indexElements( void )
    {
        int nV = 0;
        for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
        {
            v->index = nV;
            nV++;
        }

        int nE = 0;
        for( EdgeIter e = edges.begin(); e != edges.end(); e++ )
        {
            e->index = nE;
            nE++;
        }

        int nF = 0;
        for( FaceIter f = faces.begin(); f != faces.end(); f++ )
        {
            f->index = nF;
            nF++;
        }
    }

    double Mesh::area( void ) const
    {
        double sum = 0.0;
        for( FaceCIter f = faces.begin();
                f != faces.end();
                f ++ )
        {
            sum += f->area();
        }
        return sum;
    }

    double Mesh::meanDiedralAngle() const
    {
        assert( !edges.empty() );

        double accum = 0.0;
        for (const Edge& edge : edges)
        {
            if (edge.he->onBoundary || edge.he->flip->onBoundary) continue;
            const Vector aa = edge.he->face->normal.unit();
            const Vector bb = edge.he->flip->face->normal.unit();
            const double scalar_product = dot(aa,bb);
            const double angle = acos(scalar_product > 1 ? 1 : scalar_product < -1 ? -1 : scalar_product);
            assert( std::isfinite(angle) );
            accum += angle;
        }

        return accum/edges.size();
    }

    double Mesh::meanEdgeLength( void  ) const
    {
        assert( !edges.empty() );

        double sum = 0;
        for( EdgeCIter e = edges.begin();
                e != edges.end();
                e ++)
        {
            VertexIter v0 = e->he->vertex;
            VertexIter v1 = e->he->flip->vertex;
            sum += (v0->position - v1->position).norm();
        }

        return sum / edges.size();
    }
}

