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

#ifndef integrate_grid_function_hpp
#define integrate_grid_function_hpp

#include <common/armadillo_fwd.hpp>

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType> class GridFunction;
class GridSegment;

//! Return the integral of \p gridFunction over the grid on which it is defined.
template <typename BasisFunctionType, typename ResultType>
arma::Col<ResultType> integrateGridFunction(
    const GridFunction<BasisFunctionType, ResultType>& gridFunction);

//! Return the integral of \p gridFunction over the segment \p gridSegment 
//! of the grid on which it is defined.
template <typename BasisFunctionType, typename ResultType>
arma::Col<ResultType> integrateGridFunctionOnSegment(
    const GridFunction<BasisFunctionType, ResultType>& gridFunction, 
    const GridSegment &gridSegment);

} // end namespace Bempp

#endif
