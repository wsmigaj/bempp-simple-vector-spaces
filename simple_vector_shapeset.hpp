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

#ifndef simple_vector_shapeset_hpp
#define simple_vector_shapeset_hpp

#include "fiber/basis.hpp"
#include "fiber/basis_data.hpp"

#include <algorithm>

namespace Fiber
{

template <typename ValueType, int dim>
class SimpleVectorShapeset : public Basis<ValueType>
{
public:
    typedef typename Basis<ValueType>::CoordinateType CoordinateType;

    SimpleVectorShapeset(const shared_ptr<Shapeset<ValueType> > &scalarShapeset) :
        m_scalarShapeset(scalarShapeset)
    {}

    virtual int size() const {
        return m_scalarShapeset->size() * dim;
    }

    virtual int order() const {
        return m_scalarShapeset->order();
    }

    virtual void evaluate(size_t what,
                          const arma::Mat<CoordinateType>& points,
                          LocalDofIndex localDofIndex,
                          BasisData<ValueType>& data) const {
        const size_t pointCount = points.n_cols;
        // TODO: perhaps cache this as thread-local storage to avoid frequent
        // allocaltion/deallocation
        BasisData<ValueType> scalarData;
        const LocalDofIndex scalarLocalDofIndex = 
            localDofIndex == ALL_DOFS ? ALL_DOFS : localDofIndex / 3;
        m_scalarShapeset->evaluate(what, points, scalarLocalDofIndex, scalarData);
        if (localDofIndex == ALL_DOFS)
        {
            if (what & VALUES)
            {
                assert(scalarData.values.extent(0) == 1);
                assert(scalarData.values.extent(2) == pointCount);
                const size_t scalarDofCount = scalarData.values.extent(1);
                data.values.set_size(dim, scalarDofCount * dim, pointCount);
                std::fill(data.values.begin(), data.values.end(), 0);
                for (size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex)
                    for (size_t scalarDofIndex = 0; scalarDofIndex < scalarDofCount; ++scalarDofIndex)
                        for (size_t component = 0; component < dim; ++component)
                            data.values(component, scalarDofIndex * dim + component, pointIndex) = 
                                scalarData.values(0, scalarDofIndex, pointIndex);
            }
            if (what & DERIVATIVES)
            {
                assert(scalarData.derivatives.extent(0) == 1);
                assert(scalarData.derivatives.extent(3) == pointCount);
                const size_t scalarDimCount = scalarData.derivatives.extent(1);
                const size_t scalarDofCount = scalarData.derivatives.extent(2);
                data.derivatives.set_size(dim, scalarDimCount,
                                          scalarDofCount * dim, pointCount);
                std::fill(data.derivatives.begin(), data.derivatives.end(), 0);
                for (size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex)
                    for (size_t scalarDofIndex = 0; scalarDofIndex < scalarDofCount; ++scalarDofIndex)
                        for (size_t scalarDimIndex = 0; scalarDimIndex < scalarDimCount; ++scalarDimIndex)
                            for (size_t component = 0; component < dim; ++component)
                                data.derivatives(component, scalarDimIndex, 
                                            scalarDofIndex * dim + component, pointIndex) = 
                                    scalarData.derivatives(0, scalarDimIndex, 
                                                           scalarDofIndex, pointIndex);
            }
        }
        else
        {
            const size_t component = localDofIndex % dim;
            if (what & VALUES)
            {
                assert(scalarData.values.extent(0) == 1);
                assert(scalarData.values.extent(1) == 1);
                assert(scalarData.values.extent(2) == pointCount);
                data.values.set_size(dim, 1, pointCount);
                std::fill(data.values.begin(), data.values.end(), 0);
                for (size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex)
                    data.values(component, 0, pointIndex) = 
                        scalarData.values(0, 0, pointIndex);
            }
            if (what & DERIVATIVES)
            {
                assert(scalarData.derivatives.extent(0) == 1);
                assert(scalarData.derivatives.extent(2) == 1);
                assert(scalarData.derivatives.extent(3) == pointCount);
                const size_t scalarDimCount = scalarData.derivatives.extent(1);
                data.derivatives.set_size(dim, scalarDimCount, 1, pointCount);
                std::fill(data.derivatives.begin(), data.derivatives.end(), 0);
                for (size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex)
                    for (size_t scalarDimIndex = 0; scalarDimIndex < scalarDimCount; ++scalarDimIndex)
                        data.derivatives(component, scalarDimIndex, 0, pointIndex) = 
                            scalarData.derivatives(0, scalarDimIndex, 0, pointIndex);
            }
        }
    }

private:
    shared_ptr<Shapeset<ValueType> > m_scalarShapeset;
};

} // namespace Fiber

#endif
