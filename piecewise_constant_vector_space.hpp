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

#ifndef piecewise_constant_vector_space_hpp
#define piecewise_constant_vector_space_hpp

#include "simple_vector_space.hpp"

namespace Bempp
{

class GridSegment;

template <typename BasisFunctionType, int codomainDim>
class PiecewiseConstantVectorSpace : public SimpleVectorSpace<BasisFunctionType, codomainDim>
{
public:
    typedef typename Space<BasisFunctionType>::CoordinateType CoordinateType;

    enum { PIECEWISE_CONSTANT_VECTOR_BASE = 100 };

    /** \brief Constructor.
     *
     *  Construct a space of piecewise constant vector functions with \p codomainDim
     *  components defined on the grid \p grid.
     *
     *  An exception is thrown if \p grid is a null pointer.
     */
    explicit PiecewiseConstantVectorSpace(const shared_ptr<const Grid>& grid);

    /** \brief Constructor.
     *
     *  Construct a space of piecewise constant vector functions with \p
     *  codomainDim components defined on the elements of the grid \p grid
     *  belonging to the segment \p segment.
     *
     *  An exception is thrown if \p grid is a null pointer.
     */
    PiecewiseConstantVectorSpace(const shared_ptr<const Grid>& grid,
                                 const GridSegment& segment);

    virtual shared_ptr<const Space<BasisFunctionType> > discontinuousSpace(
        const shared_ptr<const Space<BasisFunctionType> >& self) const;

    virtual const Fiber::Shapeset<BasisFunctionType>& shapeset(
        const Entity<0>& element) const;

    virtual shared_ptr<const Space<BasisFunctionType> > barycentricSpace(
        const shared_ptr<const Space<BasisFunctionType> >& self) const;

    virtual SpaceIdentifier spaceIdentifier() const;

    virtual bool spaceIsCompatible(const Space<BasisFunctionType>& other) const;

private:
    /** \cond PRIVATE */
    shared_ptr<Fiber::Shapeset<BasisFunctionType> > m_shapeset;
    /** \endcond*/
};

} // namespace Bempp
#endif
