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

#ifndef piecewise_linear_continuous_vector_space_hpp
#define piecewise_linear_continuous_vector_space_hpp

#include "piecewise_linear_vector_space.hpp"

#include "grid/grid_segment.hpp"

#include <memory>
#include <tbb/mutex.h>

namespace Bempp
{

class GridSegment;

template <typename BasisFunctionType, int codomainDim>
class PiecewiseLinearContinuousVectorSpace : 
        public PiecewiseLinearVectorSpace<BasisFunctionType, codomainDim>
{
public:
    typedef typename Space<BasisFunctionType>::CoordinateType CoordinateType;

    enum { PIECEWISE_LINEAR_CONTINUOUS_VECTOR_BASE = 110 };

    /** \brief Constructor.
     *
     *  Construct a space of piecewise linear, continuous vector functions with
     *  \p codomainDim components defined on the grid \p grid.
     *
     *  An exception is thrown if \p grid is a null pointer.
     */
    explicit PiecewiseLinearContinuousVectorSpace(const shared_ptr<const Grid>& grid);

    /** \brief Constructor.
     *
     *  Construct a space of piecewise linear, continuous vector functions
     *  defined on the segment \p segment of the grid \p grid. More precisely,
     *  the space will encompass those basis functions that are associated with
     *  vertices belonging to \p segment. If \p strictlyOnSegment is \c true,
     *  the support of the basis functions is truncated to the elements that
     *  belong to \p segment, too; in this case, the space may in fact contain
     *  discontinuous basis functions when considered on the whole \p grid,
     *  although the basis functions will be continuous when considered on the
     *  chosen grid segment.
     *
     *  An exception is thrown if \p grid is a null pointer.
     */
    PiecewiseLinearContinuousVectorSpace(const shared_ptr<const Grid>& grid,
                                         const GridSegment& segment,
                                         bool strictlyOnSegment = false);

    virtual ~PiecewiseLinearContinuousVectorSpace();

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
    mutable shared_ptr<Space<BasisFunctionType> > m_discontinuousSpace;
    mutable tbb::mutex m_discontinuousSpaceMutex;
    /** \endcond */
};

} // namespace Bempp
#endif
