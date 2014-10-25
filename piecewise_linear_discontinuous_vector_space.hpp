// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef piecewise_linear_discontinuous_vector_space_hpp
#define piecewise_linear_discontinuous_vector_space_hpp

#include "piecewise_linear_vector_space.hpp"

#include "grid/grid_segment.hpp"

#include <memory>
#include <tbb/mutex.h>

namespace Bempp
{

class GridSegment;

template <typename BasisFunctionType, int codomainDim>
class PiecewiseLinearDiscontinuousVectorSpace : 
        public PiecewiseLinearVectorSpace<BasisFunctionType, codomainDim>
{
public:
    typedef typename Space<BasisFunctionType>::CoordinateType CoordinateType;

    enum { PIECEWISE_LINEAR_DISCONTINUOUS_VECTOR_BASE = 120 };

    /** \brief Constructor.
     *
     *  Construct a space of piecewise linear, discontinuous vector functions with
     *  \p codomainDim components defined on the grid \p grid.
     *
     *  An exception is thrown if \p grid is a null pointer.
     */
    explicit PiecewiseLinearDiscontinuousVectorSpace(const shared_ptr<const Grid>& grid);

    /** \brief Constructor.
     *
     *  Construct a space of piecewise linear, not necessarily continuous,
     *  vector functions defined on the segment \p segment of the grid \p grid.
     *  If \p strictlyOnSegment is set to \c false (default), the space will
     *  include all basis functions associated with vertices belonging to \p
     *  segment, regardless of whether the elements on which these functions
     *  are defined belong themselves to \p segment. In consequence, the
     *  resulting space will be (in the mathematical sense) a superset of a
     *  PiecewiseLinearContinuousVectorSpace defined on the same segment. If \p
     *  strictlyOnSegment is set to \c true, the space will only include basis
     *  functions defined on elements belonging to \p segment.
     *
     *  An exception is thrown if \p grid is a null pointer.
     */
    PiecewiseLinearDiscontinuousVectorSpace(const shared_ptr<const Grid>& grid,
                                         const GridSegment& segment,
                                         bool strictlyOnSegment = false);

    virtual shared_ptr<const Space<BasisFunctionType> > discontinuousSpace(
        const shared_ptr<const Space<BasisFunctionType> >& self) const;

    virtual shared_ptr<const Space<BasisFunctionType> > barycentricSpace(
        const shared_ptr<const Space<BasisFunctionType> >& self) const;

    virtual SpaceIdentifier spaceIdentifier() const;

    virtual bool spaceIsCompatible(const Space<BasisFunctionType>& other) const;

private:
    /** \cond PRIVATE */
    GridSegment m_segment;
    bool m_strictlyOnSegment;
    /** \endcond */
};

} // namespace Bempp
#endif
