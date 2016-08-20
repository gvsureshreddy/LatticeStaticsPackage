// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include "manually_grid_tools.h"

#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/filtered_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria_base.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/numerics/matrix_tools.h>

#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <cmath>
#include <numeric>
#include <list>
#include <set>

DEAL_II_NAMESPACE_OPEN

namespace ManualGridTools
{
   using namespace dealii;
  template<typename FaceIterator>
  inline bool
  orthogonal_equality (std::bitset<3>     &orientation,
                       const FaceIterator &face1,
                       const FaceIterator &face2,
                       const int          direction,
                       const Tensor<1,FaceIterator::AccessorType::space_dimension> &offset,
                       const FullMatrix<double> &matrix)
  {
    Assert(matrix.m() == matrix.n(),
           ExcMessage("The supplied matrix must be a square matrix"));

    static const int dim = FaceIterator::AccessorType::dimension;

    // Do a full matching of the face vertices:

    std_cxx11::
    array<unsigned int, GeometryInfo<dim>::vertices_per_face> matching;

    std::set<unsigned int> face2_vertices;
    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face; ++i)
      face2_vertices.insert(i);

    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face; ++i)
      {
        for (std::set<unsigned int>::iterator it = face2_vertices.begin();
             it != face2_vertices.end();
             it++)
          {
            if (orthogonal_equality(face1->vertex(i),face2->vertex(*it),
                                    direction, offset, matrix))
              {
                matching[i] = *it;
                face2_vertices.erase(it);
                break; // jump out of the innermost loop
              }
          }
      }

    // And finally, a lookup to determine the ordering bitmask:
    if (face2_vertices.empty())
      orientation = GridTools::OrientationLookupTable<dim>::lookup(matching);

    return face2_vertices.empty();
  }



  template<typename FaceIterator>
  inline bool
  orthogonal_equality (const FaceIterator &face1,
                       const FaceIterator &face2,
                       const int          direction,
                       const Tensor<1,FaceIterator::AccessorType::space_dimension> &offset,
                       const FullMatrix<double> &matrix)
  {
    // Call the function above with a dummy orientation array
    std::bitset<3> dummy;
    return orthogonal_equality (dummy, face1, face2, direction, offset, matrix);
  }



  /*
   * Internally used in collect_periodic_faces
   */
  template<typename CellIterator>
  void
  match_periodic_face_pairs
  (std::set<std::pair<CellIterator, unsigned int> > &pairs1,
   std::set<std::pair<typename identity<CellIterator>::type, unsigned int> > &pairs2,
   const int                                        direction,
   std::vector<GridTools::PeriodicFacePair<CellIterator> >     &matched_pairs,
   const dealii::Tensor<1,CellIterator::AccessorType::space_dimension> &offset,
   const FullMatrix<double>                         &matrix)
  {
    static const int space_dim = CellIterator::AccessorType::space_dimension;
    (void)space_dim;
    Assert (0<=direction && direction<space_dim,
            ExcIndexRange (direction, 0, space_dim));

    Assert (pairs1.size() == pairs2.size(),
            ExcMessage ("Unmatched faces on periodic boundaries"));

    unsigned int n_matches = 0;

    // Match with a complexity of O(n^2). This could be improved...
    std::bitset<3> orientation;
    typedef typename std::set
    <std::pair<CellIterator, unsigned int> >::const_iterator PairIterator;
    for (PairIterator it1 = pairs1.begin(); it1 != pairs1.end(); ++it1)
      {
        for (PairIterator it2 = pairs2.begin(); it2 != pairs2.end(); ++it2)
          {
            const CellIterator cell1 = it1->first;
            const CellIterator cell2 = it2->first;
            const unsigned int face_idx1 = it1->second;
            const unsigned int face_idx2 = it2->second;
            if (ManualGridTools::orthogonal_equality(orientation,
                                               cell1->face(face_idx1),
                                               cell2->face(face_idx2),
                                               direction, offset,
                                               matrix))
              {
                // We have a match, so insert the matching pairs and
                // remove the matched cell in pairs2 to speed up the
                // matching:
                const GridTools::PeriodicFacePair<CellIterator> matched_face =
                {
                  {cell1, cell2},
                  {face_idx1, face_idx2},
                  orientation,
                  matrix
                };
                matched_pairs.push_back(matched_face);
                pairs2.erase(it2);
                ++n_matches;
                break;
              }
          }
      }

    //Assure that all faces are matched
    AssertThrow (n_matches == pairs1.size() && pairs2.size() == 0,
                 ExcMessage ("Unmatched faces on periodic boundaries"));
  }



  template<typename MeshType>
  void
  collect_periodic_faces
  (const MeshType                        &mesh,
   const types::boundary_id               b_id1,
   const types::boundary_id               b_id2,
   const int                              direction,
   std::vector<GridTools::PeriodicFacePair<typename MeshType::cell_iterator> > &matched_pairs,
   const Tensor<1,MeshType::space_dimension> &offset,
   const FullMatrix<double>              &matrix)
  {
    static const int dim = MeshType::dimension;
    static const int space_dim = MeshType::space_dimension;
    (void)dim;
    (void)space_dim;
    Assert (0<=direction && direction<space_dim,
            ExcIndexRange (direction, 0, space_dim));

    // Loop over all cells on the highest level and collect all boundary
    // faces belonging to b_id1 and b_id2:

    std::set<std::pair<typename MeshType::cell_iterator, unsigned int> > pairs1;
    std::set<std::pair<typename MeshType::cell_iterator, unsigned int> > pairs2;

    for (typename MeshType::active_cell_iterator cell = mesh.begin_active();
         cell != mesh.end(); ++cell)
      {
        for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i)
          {
            const typename MeshType::face_iterator face = cell->face(i);
            if (face->at_boundary() && face->boundary_id() == b_id1)
              {
                const std::pair<typename MeshType::cell_iterator, unsigned int> pair1
                  = std::make_pair(cell, i);
                pairs1.insert(pair1);
              }

            if (face->at_boundary() && face->boundary_id() == b_id2)
              {
                const std::pair<typename MeshType::cell_iterator, unsigned int> pair2
                  = std::make_pair(cell, i);
                pairs2.insert(pair2);
              }
          }
      }

    Assert (pairs1.size() == pairs2.size(),
            ExcMessage ("Unmatched faces on periodic boundaries"));

    // and call match_periodic_face_pairs that does the actual matching:
    match_periodic_face_pairs(pairs1, pairs2, direction, matched_pairs, offset,
                              matrix);
  }



  template<typename MeshType>
  void
  collect_periodic_faces
  (const MeshType                        &mesh,
   const types::boundary_id               b_id,
   const int                              direction,
   std::vector<GridTools::PeriodicFacePair<typename MeshType::cell_iterator> > &matched_pairs,
   const Tensor<1,MeshType::space_dimension> &offset,
   const FullMatrix<double>              &matrix)
  {
    static const int dim = MeshType::dimension;
    static const int space_dim = MeshType::space_dimension;
    (void)dim;
    (void)space_dim;
    Assert (0<=direction && direction<space_dim,
            ExcIndexRange (direction, 0, space_dim));

    Assert(dim == space_dim,
           ExcNotImplemented());

    // Loop over all cells on the highest level and collect all boundary
    // faces 2*direction and 2*direction*1:

    std::set<std::pair<typename MeshType::cell_iterator, unsigned int> > pairs1;
    std::set<std::pair<typename MeshType::cell_iterator, unsigned int> > pairs2;

    for (typename MeshType::cell_iterator cell = mesh.begin(0);
         cell != mesh.end(0); ++cell)
      {
        const typename MeshType::face_iterator face_1 = cell->face(2*direction);
        const typename MeshType::face_iterator face_2 = cell->face(2*direction+1);

        if (face_1->at_boundary() && face_1->boundary_id() == b_id)
          {
            const std::pair<typename MeshType::cell_iterator, unsigned int> pair1
              = std::make_pair(cell, 2*direction);
            pairs1.insert(pair1);
          }

        if (face_2->at_boundary() && face_2->boundary_id() == b_id)
          {
            const std::pair<typename MeshType::cell_iterator, unsigned int> pair2
              = std::make_pair(cell, 2*direction+1);
            pairs2.insert(pair2);
          }
      }

    Assert (pairs1.size() == pairs2.size(),
            ExcMessage ("Unmatched faces on periodic boundaries"));

    // and call match_periodic_face_pairs that does the actual matching:
    match_periodic_face_pairs(pairs1, pairs2, direction, matched_pairs, offset,
                              matrix);
  }

} /* namespace ManualGridTools */