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

#ifndef dealii__manually_grid_tools_H
#define dealii__manually_grid_tools_H


#include <deal.II/grid/grid_tools.h>

#include <bitset>
#include <list>
#include <set>

/**
 * This namespace is a collection of algorithms working on triangulations,
 * such as shifting or rotating triangulations, but also finding a cell that
 * contains a given point. See the descriptions of the individual functions
 * for more information.
 *
 * @ingroup grid
 */
namespace ManualGridTools
{
    using namespace dealii;
  template<typename FaceIterator>
  bool
  orthogonal_equality (std::bitset<3>     &orientation,
                       const FaceIterator &face1,
                       const FaceIterator &face2,
                       const int          direction,
                       const Tensor<1,FaceIterator::AccessorType::space_dimension> &offset
                       = Tensor<1,FaceIterator::AccessorType::space_dimension>(),
                       const FullMatrix<double> &matrix = FullMatrix<double>());


  /**
   * Same function as above, but doesn't return the actual orientation
   */
  template<typename FaceIterator>
  bool
  orthogonal_equality (const FaceIterator &face1,
                       const FaceIterator &face2,
                       const int          direction,
                       const Tensor<2,FaceIterator::AccessorType::space_dimension> &offset
                       = Tensor<1,FaceIterator::AccessorType::space_dimension>(),
                       const FullMatrix<double> &matrix = FullMatrix<double>());


  /**
   * This function will collect periodic face pairs on the coarsest mesh level
   * of the given @p mesh (a Triangulation or DoFHandler) and add them to the
   * vector @p matched_pairs leaving the original contents intact.
   *
   * Define a 'first' boundary as all boundary faces having boundary_id @p
   * b_id1 and a 'second' boundary consisting of all faces belonging to @p
   * b_id2.
   *
   * This function tries to match all faces belonging to the first boundary
   * with faces belonging to the second boundary with the help of
   * orthogonal_equality().
   *
   * The bitset that is returned inside of PeriodicFacePair encodes the
   * _relative_ orientation of the first face with respect to the second face,
   * see the documentation of orthogonal_equality() for further details.
   *
   * The @p direction refers to the space direction in which periodicity is
   * enforced. When maching periodic faces this vector component is ignored.
   *
   * The @p offset is a vector tangential to the faces that is added to the
   * location of vertices of the 'first' boundary when attempting to match
   * them to the corresponding vertices of the 'second' boundary. This can be
   * used to implement conditions such as $u(0,y)=u(1,y+1)$.
   *
   * Optionally, a $dim\times dim$ rotation @p matrix can be specified that
   * describes how vector valued DoFs of the first face should be modified
   * prior to constraining to the DoFs of the second face. The @p matrix is
   * used in two places. First, @p matrix will be supplied to
   * orthogonal_equality() and used for matching faces: Two vertices $v_1$ and
   * $v_2$ match if $\text{matrix}\cdot v_1 + \text{offset} - v_2$ is parallel
   * to the unit vector in unit direction @p direction. (For more details see
   * DoFTools::make_periodicity_constraints(), the glossary
   * @ref GlossPeriodicConstraints "glossary entry on periodic conditions"
   * and step-45). Second, @p matrix will be stored in the PeriodicFacePair
   * collection @p matched_pairs for further use.
   *
   * @tparam MeshType A type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   *
   * @note The created std::vector can be used in
   * DoFTools::make_periodicity_constraints() and in
   * parallel::distributed::Triangulation::add_periodicity() to enforce
   * periodicity algebraically.
   *
   * @note Because elements will be added to @p matched_pairs (and existing
   * entries will be preserved), it is possible to call this function several
   * times with different boundary ids to generate a vector with all periodic
   * pairs.
   *
   * @author Daniel Arndt, Matthias Maier, 2013 - 2015
   */
  template <typename MeshType>
  void
  collect_periodic_faces
  (const MeshType                            &mesh,
   const types::boundary_id                   b_id1,
   const types::boundary_id                   b_id2,
   const int                                  direction,
   std::vector<GridTools::PeriodicFacePair<typename MeshType::cell_iterator> > &matched_pairs,
   const Tensor<1,MeshType::space_dimension> &offset = dealii::Tensor<1,MeshType::space_dimension>(),
   const FullMatrix<double>                  &matrix = FullMatrix<double>());
//   {
//    static const int dim = MeshType::dimension;
//    static const int space_dim = MeshType::space_dimension;
//    (void)dim;
//    (void)space_dim;
//    Assert (0<=direction && direction<space_dim,
//            ExcIndexRange (direction, 0, space_dim));
//
//    // Loop over all cells on the highest level and collect all boundary
//    // faces belonging to b_id1 and b_id2:
//
//    std::set<std::pair<typename MeshType::cell_iterator, unsigned int> > pairs1;
//    std::set<std::pair<typename MeshType::cell_iterator, unsigned int> > pairs2;
//
//    for (typename MeshType::active_cell_iterator cell = mesh.begin_active();
//         cell != mesh.end(); ++cell)
//      {
//        for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i)
//          {
//            const typename MeshType::face_iterator face = cell->face(i);
//            if (face->at_boundary() && face->boundary_id() == b_id1)
//              {
//                const std::pair<typename MeshType::cell_iterator, unsigned int> pair1
//                  = std::make_pair(cell, i);
//                pairs1.insert(pair1);
//              }
//
//            if (face->at_boundary() && face->boundary_id() == b_id2)
//              {
//                const std::pair<typename MeshType::cell_iterator, unsigned int> pair2
//                  = std::make_pair(cell, i);
//                pairs2.insert(pair2);
//              }
//          }
//      }
//
//    Assert (pairs1.size() == pairs2.size(),
//            ExcMessage ("Unmatched faces on periodic boundaries"));
//
//    // and call match_periodic_face_pairs that does the actual matching:
//    //match_periodic_face_pairs(pairs1, pairs2, direction, matched_pairs, offset,
//    //                          matrix);
//  }
  
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
   const FullMatrix<double>                         &matrix);
//  {
//    static const int space_dim = CellIterator::AccessorType::space_dimension;
//    (void)space_dim;
//    Assert (0<=direction && direction<space_dim,
//            ExcIndexRange (direction, 0, space_dim));
//
//    Assert (pairs1.size() == pairs2.size(),
//            ExcMessage ("Unmatched faces on periodic boundaries"));
//
//    unsigned int n_matches = 0;
//
//    // Match with a complexity of O(n^2). This could be improved...
//    std::bitset<3> orientation;
//    typedef typename std::set
//    <std::pair<CellIterator, unsigned int> >::const_iterator PairIterator;
//    for (PairIterator it1 = pairs1.begin(); it1 != pairs1.end(); ++it1)
//      {
//        for (PairIterator it2 = pairs2.begin(); it2 != pairs2.end(); ++it2)
//          {
//            const CellIterator cell1 = it1->first;
//            const CellIterator cell2 = it2->first;
//            const unsigned int face_idx1 = it1->second;
//            const unsigned int face_idx2 = it2->second;
//            if (ManualGridTools::orthogonal_equality(orientation,
//                                               cell1->face(face_idx1),
//                                               cell2->face(face_idx2),
//                                               direction, offset,
//                                               matrix))
//              {
//                // We have a match, so insert the matching pairs and
//                // remove the matched cell in pairs2 to speed up the
//                // matching:
//                const GridTools::PeriodicFacePair<CellIterator> matched_face =
//                {
//                  {cell1, cell2},
//                  {face_idx1, face_idx2},
//                  orientation,
//                  matrix
//                };
//                matched_pairs.push_back(matched_face);
//                pairs2.erase(it2);
//                ++n_matches;
//                break;
//              }
//          }
//      }
//
//    //Assure that all faces are matched
//    AssertThrow (n_matches == pairs1.size() && pairs2.size() == 0,
//                 ExcMessage ("Unmatched faces on periodic boundaries"));
//  }


  /**
   * This compatibility version of collect_periodic_face_pairs() only works on
   * grids with cells in
   * @ref GlossFaceOrientation "standard orientation".
   *
   * Instead of defining a 'first' and 'second' boundary with the help of two
   * boundary_ids this function defines a 'left' boundary as all faces with
   * local face index <code>2*dimension</code> and boundary indicator @p b_id
   * and, similarly, a 'right' boundary consisting of all face with local face
   * index <code>2*dimension+1</code> and boundary indicator @p b_id.
   *
   * This function will collect periodic face pairs on the coarsest mesh level
   * and add them to @p matched_pairs leaving the original contents intact.
   *
   * See above function for further details.
   *
   * @note This version of collect_periodic_face_pairs() will not work on
   * meshes with cells not in
   * @ref GlossFaceOrientation "standard orientation".
   *
   * @author Daniel Arndt, Matthias Maier, 2013 - 2015
   */
  template <typename MeshType>
  void
  collect_periodic_faces
  (const MeshType                                    &mesh,
   const types::boundary_id                           b_id,
   const int                                          direction,
   std::vector<GridTools::PeriodicFacePair<typename MeshType::cell_iterator> > &matched_pairs,
   const dealii::Tensor<1,MeshType::space_dimension> &offset = dealii::Tensor<1,MeshType::space_dimension>(),
   const FullMatrix<double>                          &matrix = FullMatrix<double>());

} /*namespace ManualGridTools*/
/* end of #ifndef dealii__manually_grid_tools_H */
#endif
/*----------------------------   manually_grid_tools.h     ---------------------------*/
