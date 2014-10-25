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

#ifndef simple_vector_function_value_functor_hpp
#define simple_vector_function_value_functor_hpp

#include "common/common.hpp"

#include "fiber/basis_data.hpp"
#include "fiber/geometrical_data.hpp"
#include "fiber/collection_of_3d_arrays.hpp"
#include "fiber/shape_transformation_functor_wrappers.hpp"

namespace Fiber
{

template <typename CoordinateType_, int dim>
class SimpleVectorFunctionValueElementaryFunctor
{
public:
    typedef CoordinateType_ CoordinateType;

    int argumentDimension() const { return dim; }
    int resultDimension() const { return dim; }

    void addDependencies(size_t& basisDeps, size_t& geomDeps) const {
        basisDeps |= VALUES;
    }

    template <typename ValueType>
    void evaluate(
            const ConstBasisDataSlice<ValueType>& basisData,
            const ConstGeometricalDataSlice<CoordinateType>& geomData,
            _1dSliceOf3dArray<ValueType>& result) const {
        assert(basisData.componentCount() == argumentDimension());
        assert(result.extent(0) == resultDimension());
        for (size_t i = 0; i < dim; ++i)
            result(i) = basisData.values(i);
    }
};

// Note: in C++11 we'll be able to make a "template typedef", or more precisely
// a using declaration, instead of this spurious inheritance
/** \ingroup functors
 *  \brief Functor calculating the value of a scalar basis function. */
template <typename CoordinateType_, int dim>
class SimpleVectorFunctionValueFunctor :
        public ElementaryShapeTransformationFunctorWrapper<
    SimpleVectorFunctionValueElementaryFunctor<CoordinateType_, dim> >
{
public:
    typedef CoordinateType_ CoordinateType;
};

} // namespace Fiber

#endif
