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

#include "piecewise_linear_continuous_vector_space.hpp"

#include "simple_vector_shapeset.hpp"

#include "piecewise_linear_discontinuous_vector_space.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "fiber/explicit_instantiation.hpp"

#include <boost/make_shared.hpp>

namespace Bempp
{

template <typename BasisFunctionType, int codomainDim>
PiecewiseLinearContinuousVectorSpace<BasisFunctionType, codomainDim>::PiecewiseLinearContinuousVectorSpace(
    const shared_ptr<const Grid>& grid) :
    PiecewiseLinearVectorSpace<BasisFunctionType, codomainDim>( 
        boost::make_shared<PiecewiseLinearContinuousScalarSpace<BasisFunctionType> >(grid)),
    m_segment(GridSegment::wholeGrid(*grid)),
    m_strictlyOnSegment(false)
{
}

template <typename BasisFunctionType, int codomainDim>
PiecewiseLinearContinuousVectorSpace<BasisFunctionType, codomainDim>::PiecewiseLinearContinuousVectorSpace(
    const shared_ptr<const Grid>& grid,
    const GridSegment& segment,
    bool strictlyOnSegment) :
    PiecewiseLinearVectorSpace<BasisFunctionType, codomainDim>( 
        boost::make_shared<PiecewiseLinearContinuousScalarSpace<BasisFunctionType> >(
            grid, segment, strictlyOnSegment)),
    m_segment(segment),
    m_strictlyOnSegment(strictlyOnSegment)
{
}

template <typename BasisFunctionType, int codomainDim>
PiecewiseLinearContinuousVectorSpace<BasisFunctionType, codomainDim>::~PiecewiseLinearContinuousVectorSpace()
{
}

template <typename BasisFunctionType, int codomainDim>
shared_ptr<const Space<BasisFunctionType> >
PiecewiseLinearContinuousVectorSpace<BasisFunctionType, codomainDim>::discontinuousSpace(
    const shared_ptr<const Space<BasisFunctionType> >& self) const
{
    if (!m_discontinuousSpace) {
        tbb::mutex::scoped_lock lock(m_discontinuousSpaceMutex);
        typedef PiecewiseLinearDiscontinuousVectorSpace<BasisFunctionType, codomainDim>
                DiscontinuousSpace;
        if (!m_discontinuousSpace)
            m_discontinuousSpace.reset(
                        new DiscontinuousSpace(this->grid(), m_segment,
                                               m_strictlyOnSegment));
    }
    return m_discontinuousSpace;
}

template <typename BasisFunctionType, int codomainDim>
shared_ptr<const Space<BasisFunctionType> > 
PiecewiseLinearContinuousVectorSpace<BasisFunctionType, codomainDim>::barycentricSpace(
    const shared_ptr<const Space<BasisFunctionType> >& self) const
{
    throw std::runtime_error("PiecewiseLinearContinuousVectorSpace::barycentricSpace(): "
                             "not implemented yet");
}

template <typename BasisFunctionType, int codomainDim>
SpaceIdentifier
PiecewiseLinearContinuousVectorSpace<BasisFunctionType, codomainDim>::spaceIdentifier() const
{
    // Ugly hack. This function, if it needs to exist at all, should return
    // an int, not a SpaceIdentifier
    return static_cast<SpaceIdentifier>(
        PIECEWISE_LINEAR_CONTINUOUS_VECTOR_BASE + codomainDim);
}

template <typename BasisFunctionType, int codomainDim>
bool
PiecewiseLinearContinuousVectorSpace<BasisFunctionType, codomainDim>::spaceIsCompatible(
    const Space<BasisFunctionType>& other) const
{
    return other.grid().get() == this->grid().get() && 
        other.spaceIdentifier() == this->spaceIdentifier();
}

#define INSTANTIATE_PIECEWISE_LINEAR_CONTINUOUS_VECTOR_SPACE(BASIS) \
    template class PiecewiseLinearContinuousVectorSpace< BASIS, 3 >;
FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_PIECEWISE_LINEAR_CONTINUOUS_VECTOR_SPACE);

} // namespace Bempp

