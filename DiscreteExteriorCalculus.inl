#include "DiscreteExteriorCalculus.h"

namespace DDG
{
   template<typename T> T clamp(const T& t)
   {
       static const T m = 1e6;
       assert( !std::isnan(t) );
       return t > m ? m : t < -m ? -m : t;
   }

   template< class T >
   void Diag<T>::build( const DenseMatrix<T>& vector, SparseMatrix<T>& matrix )
   {
       matrix = SparseMatrix<Real>(vector.length(), vector.length());

       for (int index=0; index<vector.length(); ++index)
           matrix(index, index) = vector(index);
   }

   template <class T>
   void Identity0Form<T> :: build( const Mesh& mesh,
                                    SparseMatrix<T>& id0 )
   {
      int nV = mesh.vertices.size();

      id0 = SparseMatrix<T>( nV, nV );

      for( VertexCIter v  = mesh.vertices.begin();
                       v != mesh.vertices.end();
                       v ++ )
      {
         int i = v->index;
         id0( i, i ) = 1;
      }
   }

   template <class T>
   void Identity1Form<T> :: build( const Mesh& mesh,
                                    SparseMatrix<T>& id1 )
   {
      int nE = mesh.edges.size();

      id1 = SparseMatrix<T>( nE, nE );

      for( EdgeCIter e  = mesh.edges.begin();
                       e != mesh.edges.end();
                       e ++ )
      {
         int i = e->index;
         id1( i, i ) = 1;
      }
   }

   template <class T>
   void Identity2Form<T> :: build( const Mesh& mesh,
                                    SparseMatrix<T>& id2 )
   {
      int nF = mesh.faces.size();

      id2 = SparseMatrix<T>( nF, nF );

      for( FaceCIter f  = mesh.faces.begin();
                       f != mesh.faces.end();
                       f ++ )
      {
         int i = f->index;
         id2( i, i ) = 1;
      }
   }

   template <class T>
   void HodgeStar0Form<T> :: build( const Mesh& mesh, const bool inv,
                                    SparseMatrix<T>& star0 )
   // builds a diagonal matrix mapping primal discrete 0-forms
   // to dual discrete 2-forms
   {
      int nV = mesh.vertices.size();

      star0 = SparseMatrix<T>( nV, nV );

      for( VertexCIter v  = mesh.vertices.begin();
                       v != mesh.vertices.end();
                       v ++ )
      {
         int i = v->index;
         star0( i, i ) = clamp(inv ? 1./v->area() : v->area());
      }
   }

   template <class T>
   void HodgeStar1Form<T> :: build( const Mesh& mesh, const bool inv,
                                    SparseMatrix<T>& star1 )
   // builds a diagonal matrix mapping primal discrete 1-forms
   // to dual discrete 1-forms
   {
      int nE = mesh.edges.size();

      star1 = SparseMatrix<T>( nE, nE );

      for( EdgeCIter e  = mesh.edges.begin();
                     e != mesh.edges.end();
                     e ++ )
      {
         // get the cotangents of the two angles opposite this edge
         double cotAlpha = e->he->cotan();
         double cotBeta  = e->he->flip->cotan();

         int i = e->index;
         star1( i, i ) = clamp(inv ? 2./(cotAlpha + cotBeta) : ( cotAlpha + cotBeta )/2.);
      }
   }

   template <class T>
   void HodgeStar2Form<T> :: build( const Mesh& mesh, const bool inv,
                                    SparseMatrix<T>& star2 )
   // builds a diagonal matrix mapping primal discrete 2-forms
   // to dual discrete 2-forms
   {
      int nF = mesh.faces.size();

      star2 = SparseMatrix<T>( nF, nF );

      for( FaceCIter f  = mesh.faces.begin();
                     f != mesh.faces.end();
                     f ++ )
      {
         int i = f->index;
         star2( i, i ) = clamp(inv ? f->area() : 1./f->area());
      }
   }

   template< class T >
   void ExteriorDerivative0Form<T> :: build( const Mesh& mesh,
                                             SparseMatrix<T>& d0 )
   {
      int nV = mesh.vertices.size();
      int nE = mesh.edges.size();

      d0 = SparseMatrix<T>( nE, nV );

      for( EdgeCIter e  = mesh.edges.begin();
                     e != mesh.edges.end();
                     e ++ )
      {
         // the row index is the index of the edge
         int r = e->index;

         // the column indices are the indices of the two
         // edge vertices -- orientation is determined by
         // the orientation of the edge's first half edge
         int ci = e->he->vertex->index;
         int cj = e->he->flip->vertex->index;

         d0( r, ci ) = -1.;
         d0( r, cj ) =  1.;
      }
   }

   template< class T >
   void Wedge0FormTo1Form<T> :: build( const Mesh& mesh,
                                             SparseMatrix<T>& mm )
   {
      int nV = mesh.vertices.size();
      int nE = mesh.edges.size();

      mm = SparseMatrix<T>( nE, nV );

      for( EdgeCIter e  = mesh.edges.begin();
                     e != mesh.edges.end();
                     e ++ )
      {
         // the row index is the index of the edge
         int r = e->index;

         // the column indices are the indices of the two
         // edge vertices -- orientation is determined by
         // the orientation of the edge's first half edge
         int ci = e->he->vertex->index;
         int cj = e->he->flip->vertex->index;

         mm( r, ci ) = .5;
         mm( r, cj ) = .5;
      }
   }

   template< class T >
   void ExteriorDerivative1Form<T> :: build( const Mesh& mesh,
                                             SparseMatrix<T>& d1 )
   {
      int nE = mesh.edges.size();
      int nF = mesh.faces.size();

      d1 = SparseMatrix<T>( nF, nE );

      // visit each face
      for( FaceCIter f  = mesh.faces.begin();
                     f != mesh.faces.end();
                     f ++ )
      {
         // the row index is the index of the face
         int r = f->index;

         // visit all edges of this face
         HalfEdgeCIter he = f->he;
         do
         {
            // the column index is the index of the current edge
            int c = he->edge->index;

            // relative orientation is determined by checking if
            // the current half edge is the first half edge of its
            // corresponding edge
            double s = ( he->edge->he == he ? 1. : -1. );

            // set the entry for this edge
            d1( r, c ) = s;

            he = he->next;
         }
         while( he != f->he );
      }
   }
}

