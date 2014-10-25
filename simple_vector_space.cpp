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

#include "simple_vector_space.hpp"

#include "simple_vector_function_value_functor.hpp"

#include "common/acc.hpp"

#include "fiber/explicit_instantiation.hpp"
#include "fiber/default_collection_of_shapeset_transformations.hpp"

namespace Bempp
{

/** \cond PRIVATE */
template <typename BasisFunctionType, int codomainDim>
struct SimpleVectorSpace<BasisFunctionType, codomainDim>::Impl
{
    typedef Fiber::SimpleVectorFunctionValueFunctor<CoordinateType, codomainDim>
    TransformationFunctor;

    Impl() : transformations(TransformationFunctor())
    {}

    Fiber::DefaultCollectionOfShapesetTransformations<TransformationFunctor>
    transformations;
};
/** \endcond */

template <typename BasisFunctionType, int codomainDim>
SimpleVectorSpace<BasisFunctionType, codomainDim>::SimpleVectorSpace(
    const shared_ptr<Space<BasisFunctionType> >& scalarSpace) :
    Base(scalarSpace->grid()), m_scalarSpace(scalarSpace), m_impl(new Impl)
{
    if (scalarSpace->codomainDimension() != 1)
    {
        throw std::invalid_argument("SimpleVectorSpace::SimpleVectorSpace(): "
                                    "argument must be a scalar space");
    }
}

template <typename BasisFunctionType, int codomainDim>
SimpleVectorSpace<BasisFunctionType, codomainDim>::SimpleVectorSpace(
    const SimpleVectorSpace& other) :
    Base(other), m_scalarSpace(other.m_scalarSpace), m_impl(new Impl(*other.m_impl))
{
}

template <typename BasisFunctionType, int codomainDim>
SimpleVectorSpace<BasisFunctionType, codomainDim>::~SimpleVectorSpace()
{
}

template <typename BasisFunctionType, int codomainDim>
SimpleVectorSpace<BasisFunctionType, codomainDim>&
SimpleVectorSpace<BasisFunctionType, codomainDim>::operator=(const SimpleVectorSpace& rhs)
{
    if (this != &rhs) {
        Base::operator=(rhs);
        m_scalarSpace = rhs.m_scalarSpace;
        m_impl.reset(new Impl(*rhs.m_impl));
    }
    return *this;
}

template <typename BasisFunctionType, int codomainDim>
bool SimpleVectorSpace<BasisFunctionType, codomainDim>::isDiscontinuous() const
{
    return m_scalarSpace->isDiscontinuous();
}

template <typename BasisFunctionType, int codomainDim>
bool SimpleVectorSpace<BasisFunctionType, codomainDim>::isBarycentric() const
{
    return m_scalarSpace->isBarycentric();
}

template <typename BasisFunctionType, int codomainDim>
int SimpleVectorSpace<BasisFunctionType, codomainDim>::domainDimension() const 
{
    return m_scalarSpace->domainDimension();
}

template <typename BasisFunctionType, int codomainDim>
int SimpleVectorSpace<BasisFunctionType, codomainDim>::codomainDimension() const
{
    return codomainDim;
}

template <typename BasisFunctionType, int codomainDim>
const typename SimpleVectorSpace<BasisFunctionType, codomainDim>::CollectionOfShapesetTransformations&
SimpleVectorSpace<BasisFunctionType, codomainDim>::basisFunctionValue() const 
{
    return m_impl->transformations;
}

template <typename BasisFunctionType, int codomainDim>
void SimpleVectorSpace<BasisFunctionType, codomainDim>::setElementVariant(
    const Entity<0>& element, ElementVariant variant)
{
    m_scalarSpace->setElementVariant(element, variant);
}

template <typename BasisFunctionType, int codomainDim>
ElementVariant SimpleVectorSpace<BasisFunctionType, codomainDim>::elementVariant(
    const Entity<0>& element) const
{
    return m_scalarSpace->elementVariant(element);
}

template <typename BasisFunctionType, int codomainDim>
size_t SimpleVectorSpace<BasisFunctionType, codomainDim>::flatLocalDofCount() const
{
    return m_scalarSpace->flatLocalDofCount() * codomainDim;
}

template <typename BasisFunctionType, int codomainDim>
size_t SimpleVectorSpace<BasisFunctionType, codomainDim>::globalDofCount() const
{
    return m_scalarSpace->globalDofCount() * codomainDim;
}

template <typename BasisFunctionType, int codomainDim>
void SimpleVectorSpace<BasisFunctionType, codomainDim>::getGlobalDofs(
    const Entity<0>& element,
    std::vector<GlobalDofIndex>& dofs,
    std::vector<BasisFunctionType>& localDofWeights) const
{
    m_scalarSpace->getGlobalDofs(element, dofs, localDofWeights);
    const size_t scalarDofCount = dofs.size();
    dofs.resize(codomainDim * scalarDofCount);
    for (ptrdiff_t i = scalarDofCount - 1; i >= 0; --i)
    {
        for (ptrdiff_t d = codomainDim - 1; d >= 0; --d)
        {
            acc(dofs, i * codomainDim + d) = acc(dofs, i) * codomainDim + d;
        }
    }
    localDofWeights.resize(codomainDim * scalarDofCount);
    for (ptrdiff_t i = scalarDofCount - 1; i >= 0; --i)
    {
        for (ptrdiff_t d = codomainDim - 1; d >= 0; --d)
        {
            acc(localDofWeights, i * codomainDim + d) = acc(localDofWeights, i);
        }
    }
}

template <typename BasisFunctionType, int codomainDim>
void SimpleVectorSpace<BasisFunctionType, codomainDim>::global2localDofs(
    const std::vector<GlobalDofIndex>& globalDofs,
    std::vector<std::vector<LocalDof> >& localDofs,
    std::vector<std::vector<BasisFunctionType> >& localDofWeights) const
{
    const size_t globalDofCount = globalDofs.size();
    // For efficiency (to avoid frequent reallocations), perhaps we
    // should store this vector as a thread-local member variable?
    std::vector<GlobalDofIndex> scalarGlobalDofs(globalDofCount);
    for (size_t i = 0; i < globalDofCount; ++i)
        acc(scalarGlobalDofs, i) = acc(globalDofs, i) / codomainDim;
    m_scalarSpace->global2localDofs(scalarGlobalDofs, localDofs, localDofWeights);
    for (size_t i = 0; i < globalDofCount; ++i)
    {
        const size_t component = acc(globalDofs, i) % codomainDim;
        for (size_t j = 0; j < acc(localDofs, i).size(); ++j)
            acc(acc(localDofs, i), j).dofIndex = 
                acc(acc(localDofs, i), j).dofIndex * codomainDim + component;
    }
}

template <typename BasisFunctionType, int codomainDim>
void SimpleVectorSpace<BasisFunctionType, codomainDim>::flatLocal2localDofs(
    const std::vector<FlatLocalDofIndex>& flatLocalDofs,
    std::vector<LocalDof>& localDofs) const
{
    const size_t flatLocalDofCount = flatLocalDofs.size();
    // For efficiency (to avoid frequent reallocations), perhaps we
    // should store this vector as a thread-local member variable?
    std::vector<FlatLocalDofIndex> scalarFlatLocalDofs(flatLocalDofCount);
    for (size_t i = 0; i < flatLocalDofCount; ++i)
        acc(scalarFlatLocalDofs, i) = acc(flatLocalDofs, i) / codomainDim;
    m_scalarSpace->flatLocal2localDofs(scalarFlatLocalDofs, localDofs);
    for (size_t i = 0; i < flatLocalDofCount; ++i)
    {
        const size_t component = acc(flatLocalDofs, i) % codomainDim;
        acc(localDofs, i).dofIndex =
            acc(localDofs, i).dofIndex * codomainDim + component;
    }
}

template <typename BasisFunctionType, int codomainDim>
void SimpleVectorSpace<BasisFunctionType, codomainDim>::getGlobalDofInterpolationPoints(
    arma::Mat<CoordinateType>& points) const 
{
    arma::Mat<CoordinateType> scalarSpacePoints;
    m_scalarSpace->getGlobalDofInterpolationPoints(scalarSpacePoints);
    points.set_size(scalarSpacePoints.n_rows, scalarSpacePoints.n_cols * codomainDim);
    for (size_t pointIndex = 0; pointIndex < scalarSpacePoints.n_cols; ++pointIndex)
        for (size_t component = 0; component < codomainDim; ++component)
            for (size_t dim = 0; dim < scalarSpacePoints.n_rows; ++dim)
                points(dim, pointIndex * codomainDim + component) = 
                    scalarSpacePoints(dim, pointIndex);
}

template <typename BasisFunctionType, int codomainDim>
void 
SimpleVectorSpace<BasisFunctionType, codomainDim>::getNormalsAtGlobalDofInterpolationPoints(
    arma::Mat<CoordinateType>& normals) const 
{
    arma::Mat<CoordinateType> scalarSpaceNormals;
    m_scalarSpace->getNormalsAtGlobalDofInterpolationPoints(scalarSpaceNormals);
    normals.set_size(scalarSpaceNormals.n_rows, scalarSpaceNormals.n_cols * codomainDim);
    for (size_t pointIndex = 0; pointIndex < scalarSpaceNormals.n_cols; ++pointIndex)
        for (size_t component = 0; component < codomainDim; ++component)
            for (size_t dim = 0; dim < scalarSpaceNormals.n_rows; ++dim)
                normals(dim, pointIndex * codomainDim + component) = 
                    scalarSpaceNormals(dim, pointIndex);
}

template <typename BasisFunctionType, int codomainDim>
void
SimpleVectorSpace<BasisFunctionType, codomainDim>::getGlobalDofInterpolationDirections(
    arma::Mat<CoordinateType>& directions) const 
{
    const size_t scalarDofCount = m_scalarSpace->globalDofCount();
    directions.set_size(codomainDim, codomainDim * scalarDofCount);
    directions.fill(0);
    for (size_t dofIndex = 0; dofIndex < scalarDofCount; ++dofIndex)
        for (size_t dim = 0; dim < codomainDim; ++dim)
            directions(dim, dofIndex * codomainDim + dim) = 1;
}

template <typename BasisFunctionType, int codomainDim>
void SimpleVectorSpace<BasisFunctionType, codomainDim>::getGlobalDofBoundingBoxes(
    std::vector<BoundingBox<CoordinateType> >& boundingBoxes) const 
{
    std::vector<BoundingBox<CoordinateType> > scalarDofBoundingBoxes;
    m_scalarSpace->getGlobalDofBoundingBoxes(scalarDofBoundingBoxes);
    const size_t scalarDofCount = scalarDofBoundingBoxes.size();
    boundingBoxes.resize(scalarDofCount * codomainDim);
    for (size_t dof = 0; dof < scalarDofCount; ++dof)
        for (size_t component = 0; component < codomainDim; ++component)
            acc(boundingBoxes, dof * codomainDim + component) = acc(scalarDofBoundingBoxes, dof);
}

template <typename BasisFunctionType, int codomainDim>
void SimpleVectorSpace<BasisFunctionType, codomainDim>::getFlatLocalDofBoundingBoxes(
    std::vector<BoundingBox<CoordinateType> >& boundingBoxes) const 
{
    std::vector<BoundingBox<CoordinateType> > scalarDofBoundingBoxes;
    m_scalarSpace->getFlatLocalDofBoundingBoxes(scalarDofBoundingBoxes);
    const size_t scalarDofCount = scalarDofBoundingBoxes.size();
    boundingBoxes.resize(scalarDofCount * codomainDim);
    for (size_t dof = 0; dof < scalarDofCount; ++dof)
        for (size_t component = 0; component < codomainDim; ++component)
            acc(boundingBoxes, dof * codomainDim + component) = acc(scalarDofBoundingBoxes, dof);
}

template <typename BasisFunctionType, int codomainDim>
void SimpleVectorSpace<BasisFunctionType, codomainDim>::getGlobalDofPositions(
    std::vector<Point3D<CoordinateType> >& positions) const
{
    std::vector<Point3D<CoordinateType> > scalarDofPositions;
    m_scalarSpace->getGlobalDofPositions(scalarDofPositions);
    const size_t scalarDofCount = scalarDofPositions.size();
    positions.resize(scalarDofCount * codomainDim);
    for (size_t dof = 0; dof < scalarDofCount; ++dof)
        for (size_t component = 0; component < codomainDim; ++component)
            acc(positions, dof * codomainDim + component) = acc(scalarDofPositions, dof);
}

template <typename BasisFunctionType, int codomainDim>
void SimpleVectorSpace<BasisFunctionType, codomainDim>::getFlatLocalDofPositions(
    std::vector<Point3D<CoordinateType> >& positions) const 
{
    std::vector<Point3D<CoordinateType> > scalarDofPositions;
    m_scalarSpace->getFlatLocalDofPositions(scalarDofPositions);
    const size_t scalarDofCount = scalarDofPositions.size();
    positions.resize(scalarDofCount * codomainDim);
    for (size_t dof = 0; dof < scalarDofCount; ++dof)
        for (size_t component = 0; component < codomainDim; ++component)
            acc(positions, dof * codomainDim + component) = acc(scalarDofPositions, dof);
}

template <typename BasisFunctionType, int codomainDim>
void SimpleVectorSpace<BasisFunctionType, codomainDim>::getGlobalDofNormals(
    std::vector<Point3D<CoordinateType> >& normals) const 
{
    std::vector<Point3D<CoordinateType> > scalarDofNormals;
    m_scalarSpace->getGlobalDofNormals(scalarDofNormals);
    const size_t scalarDofCount = scalarDofNormals.size();
    normals.resize(scalarDofCount * codomainDim);
    for (size_t dof = 0; dof < scalarDofCount; ++dof)
        for (size_t component = 0; component < codomainDim; ++component)
            acc(normals, dof * codomainDim + component) = acc(scalarDofNormals, dof);
}

template <typename BasisFunctionType, int codomainDim>
void SimpleVectorSpace<BasisFunctionType, codomainDim>::getFlatLocalDofNormals(
    std::vector<Point3D<CoordinateType> >& normals) const 
{
    std::vector<Point3D<CoordinateType> > scalarDofNormals;
    m_scalarSpace->getFlatLocalDofNormals(scalarDofNormals);
    const size_t scalarDofCount = scalarDofNormals.size();
    normals.resize(scalarDofCount * codomainDim);
    for (size_t dof = 0; dof < scalarDofCount; ++dof)
        for (size_t component = 0; component < codomainDim; ++component)
            acc(normals, dof * codomainDim + component) = acc(scalarDofNormals, dof);
}

template <typename BasisFunctionType, int codomainDim>
void SimpleVectorSpace<BasisFunctionType, codomainDim>::dumpClusterIds(
    const char* fileName,
    const std::vector<unsigned int>& clusterIdsOfGlobalDofs) const
{
    throw std::runtime_error("SimpleVectorSpace::dumpClusterIds(): not implemented");
}

#define INSTANTIATE_SIMPLE_VECTOR_SPACE(BASIS) \
    template class SimpleVectorSpace< BASIS, 3 >;
FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_SIMPLE_VECTOR_SPACE);

} // namespace Bempp
