#include "Face.h"
#include "Mesh.h"
#include "Vector.h"

namespace DDG
{
   double Face::area( void ) const
   {
      Vector p0 = he->vertex->position;
      Vector p1 = he->next->vertex->position;
      Vector p2 = he->next->next->vertex->position;
      return cross( p1-p0, p2-p0 ).norm() / 2.;
   }

   Vector Face::from_vertices_normal() const
   {
       const Vector n0 = he->vertex->normal;
       const Vector n1 = he->next->vertex->normal;
       const Vector n2 = he->next->next->vertex->normal;

       if (n0.norm() == 0) return Vector();
       if (n1.norm() == 0) return Vector();
       if (n2.norm() == 0) return Vector();

       return (n0+n1+n2).unit();
   }

   Vector Face::geometric_normal( void ) const
   {
      Vector p0 = he->vertex->position;
      Vector p1 = he->next->vertex->position;
      Vector p2 = he->next->next->vertex->position;

      return cross( p1-p0, p2-p0 ).unit();
   }

   bool Face::isBoundary( void ) const
   {
      return he->onBoundary;
   }

   Vector Face :: circumcenter( void ) const
   {
      Vector p0 = he->vertex->position;
      Vector p1 = he->next->vertex->position;
      Vector n = he->rotatedEdge();
      double h  = 0.5*he->cotan();
      return 0.5*(p0+p1) + h*n;
   }

   Vector Face :: barycenter( void ) const
   {
      Vector p0 = he->vertex->position;
      Vector p1 = he->next->vertex->position;
      Vector p2 = he->next->next->vertex->position;
      return (p0 + p1 + p2)/3.;
   }
}

