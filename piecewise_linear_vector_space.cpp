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

#include "piecewise_linear_vector_space.hpp"

#include "simple_vector_shapeset.hpp"

#include "space/piecewise_linear_scalar_space.hpp"
#include "fiber/linear_scalar_shapeset.hpp"
#include "fiber/explicit_instantiation.hpp"

#include <boost/make_shared.hpp>

namespace Bempp
{

template <typename BasisFunctionType, int codomainDim>
PiecewiseLinearVectorSpace<BasisFunctionType, codomainDim>::PiecewiseLinearVectorSpace(
    const shared_ptr<Space<BasisFunctionType> >& scalarSpace) :
    SimpleVectorSpace<BasisFunctionType, codomainDim>(scalarSpace),
    m_lineShapeset(
        boost::make_shared<Fiber::SimpleVectorShapeset<BasisFunctionType, codomainDim> >(
            boost::make_shared<Fiber::LinearScalarShapeset<2, BasisFunctionType> >())),
    m_triangleShapeset(
        boost::make_shared<Fiber::SimpleVectorShapeset<BasisFunctionType, codomainDim> >(
            boost::make_shared<Fiber::LinearScalarShapeset<3, BasisFunctionType> >())),
    m_quadrilateralShapeset(
        boost::make_shared<Fiber::SimpleVectorShapeset<BasisFunctionType, codomainDim> >(
            boost::make_shared<Fiber::LinearScalarShapeset<4, BasisFunctionType> >()))
{
}

template <typename BasisFunctionType, int codomainDim>
const Fiber::Shapeset<BasisFunctionType>& 
PiecewiseLinearVectorSpace<BasisFunctionType, codomainDim>::shapeset(
    const Entity<0>& element) const
{
    switch (this->elementVariant(element))
    {
    case 3:
        return *m_triangleShapeset;
    case 4:
        return *m_quadrilateralShapeset;
    case 2:
        return *m_lineShapeset;
    default:
        throw std::logic_error("PiecewiseLinearVectorSpace::shapeset(): "
                               "invalid element variant, this shouldn't happen!");
    }
}

#define INSTANTIATE_PIECEWISE_LINEAR_VECTOR_SPACE(BASIS) \
    template class PiecewiseLinearVectorSpace< BASIS, 3 >;
FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_PIECEWISE_LINEAR_VECTOR_SPACE);

} // namespace Bempp

