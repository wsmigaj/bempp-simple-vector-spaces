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

#ifndef simple_vector_space_hpp
#define simple_vector_space_hpp

#include "space/space.hpp"

#include <boost/scoped_ptr.hpp>

namespace Bempp
{

template <typename BasisFunctionType, int codomainDim>
class SimpleVectorSpace : public Space<BasisFunctionType>
{
    typedef Space<BasisFunctionType> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::CollectionOfShapesetTransformations
    CollectionOfShapesetTransformations;

    explicit SimpleVectorSpace(const shared_ptr<Space<BasisFunctionType> > &scalarSpace);

    SimpleVectorSpace(const SimpleVectorSpace& other);

    virtual ~SimpleVectorSpace();

    SimpleVectorSpace& operator=(const SimpleVectorSpace& rhs);

    virtual bool isDiscontinuous() const;

    virtual bool isBarycentric() const;

    virtual int domainDimension() const;

    virtual int codomainDimension() const;

    virtual const CollectionOfShapesetTransformations& basisFunctionValue() const;

    virtual void setElementVariant(const Entity<0>& element, ElementVariant variant);

    virtual ElementVariant elementVariant(const Entity<0>& element) const;

    virtual size_t flatLocalDofCount() const;

    virtual size_t globalDofCount() const;

    virtual void getGlobalDofs(const Entity<0>& element,
                               std::vector<GlobalDofIndex>& dofs,
                               std::vector<BasisFunctionType>& localDofWeights) const;

    virtual void global2localDofs(
            const std::vector<GlobalDofIndex>& globalDofs,
            std::vector<std::vector<LocalDof> >& localDofs,
            std::vector<std::vector<BasisFunctionType> >& localDofWeights) const;

    virtual void flatLocal2localDofs(
            const std::vector<FlatLocalDofIndex>& flatLocalDofs,
            std::vector<LocalDof>& localDofs) const;

    virtual void getGlobalDofInterpolationPoints(
        arma::Mat<CoordinateType>& points) const;

    virtual void getNormalsAtGlobalDofInterpolationPoints(
        arma::Mat<CoordinateType>& normals) const;

    virtual void getGlobalDofInterpolationDirections(
        arma::Mat<CoordinateType>& directions) const;

    virtual void getGlobalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType> >& boundingBoxes) const;  

    virtual void getFlatLocalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType> >& boundingBoxes) const;

    virtual void getGlobalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const;

    virtual void getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const;

    virtual void getGlobalDofNormals(
        std::vector<Point3D<CoordinateType> >& normals) const;

    virtual void getFlatLocalDofNormals(
        std::vector<Point3D<CoordinateType> >& normals) const;

    BEMPP_DEPRECATED virtual void dumpClusterIds(
            const char* fileName,
            const std::vector<unsigned int>& clusterIdsOfGlobalDofs) const;

protected:
    shared_ptr<Space<BasisFunctionType> > scalarSpace() { return m_scalarSpace; }

    shared_ptr<const Space<BasisFunctionType> > scalarSpace() const { return m_scalarSpace; }

private:
    /** \cond PRIVATE*/
    shared_ptr<Space<BasisFunctionType> > m_scalarSpace;

    struct Impl;
    boost::scoped_ptr<Impl> m_impl;
    /** \endcond */
};

} // namespace Bempp

#endif
