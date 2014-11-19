#include "integrate_grid_function.hpp"

#include "common/scalar_traits.hpp"
#include "assembly/grid_function.hpp"
#include "fiber/basis_data.hpp"
#include "fiber/default_single_quadrature_rule_family.hpp"
#include "fiber/explicit_instantiation.hpp"
#include "fiber/geometrical_data.hpp"
#include "fiber/shapeset.hpp"
#include "grid/entity.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/geometry.hpp"
#include "grid/grid.hpp"
#include "grid/grid_segment.hpp"
#include "grid/grid_view.hpp"
#include "space/space.hpp"

#include <iostream>
#include <map>
#include <vector>

namespace Bempp
{

namespace
{
    template <typename BasisFunctionType>
    struct ShapesetData
    {
        typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

        int functionCount;
        arma::Mat<CoordinateType> quadPoints;
        std::vector<CoordinateType> quadWeights;
        Fiber::BasisData<BasisFunctionType> basisData;
    };
}

template <typename BasisFunctionType, typename ResultType>
arma::Col<ResultType> integrateGridFunction(
    const GridFunction<BasisFunctionType, ResultType>& gridFunction)
{
    return integrateGridFunctionOnSegment(
        gridFunction, GridSegment::wholeGrid(*gridFunction.space()->grid()));
}

template <typename BasisFunctionType, typename ResultType>
arma::Col<ResultType> integrateGridFunctionOnSegment(
    const GridFunction<BasisFunctionType, ResultType>& gridFunction, 
    const GridSegment &gridSegment)
{
    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    const arma::Col<ResultType>& coeffs = gridFunction.coefficients();
    const Space<BasisFunctionType>& space = *gridFunction.space();

    const int codomainDim = space.codomainDimension();
    arma::Col<ResultType> integral(codomainDim);
    integral.fill(0.);

    const Grid& grid = *space.grid();
    std::auto_ptr<GridView> view = grid.leafView();
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    const IndexSet& indexSet = view->indexSet();

    std::vector<GlobalDofIndex> globalDofs;
    std::vector<BasisFunctionType> localDofWeights;
    std::vector<ResultType> localCoeffs;

    typedef std::map<const Fiber::Shapeset<BasisFunctionType>*,
        ShapesetData<BasisFunctionType> > ShapesetCache;
    ShapesetCache shapesetCache;

    Fiber::DefaultSingleQuadratureRuleFamily<CoordinateType> quadRuleFamily;

    Fiber::GeometricalData<CoordinateType> geomData;

    for (; !it->finished(); it->next()) {
        const Entity<0>& element = it->entity();
        if (!gridSegment.contains(0 /*codim*/, indexSet.entityIndex(element)))
            continue;

        space.getGlobalDofs(element, globalDofs, localDofWeights);
        localCoeffs.resize(globalDofs.size());

        bool anyDofsUsed = false;
        for (size_t i = 0; i < globalDofs.size(); ++i)
            if (globalDofs[i] >= 0) {
                anyDofsUsed = true;
                localCoeffs[i] = coeffs(globalDofs[i]) * localDofWeights[i];
            } else {
                localCoeffs[i] = 0.;
            }
        if (!anyDofsUsed)
            continue;

        const Geometry &geometry = element.geometry();

        const Fiber::Shapeset<BasisFunctionType>& shapeset = space.shapeset(element);
        typename ShapesetCache::const_iterator shapesetIt = 
            shapesetCache.find(&shapeset);
        if (shapesetIt == shapesetCache.end())
        {
            ShapesetData<BasisFunctionType> data;
            data.functionCount = shapeset.size();

            Fiber::SingleQuadratureDescriptor desc;
            desc.vertexCount = geometry.cornerCount();
            desc.order = shapeset.order();

            quadRuleFamily.fillQuadraturePointsAndWeights(
                desc, data.quadPoints, data.quadWeights);

            // These would need to be set differently if arbitrary functionals
            // were allowed.
            const size_t basisDataType = Fiber::VALUES;
            shapeset.evaluate(basisDataType, data.quadPoints, Fiber::ALL_DOFS, data.basisData);

            shapesetIt = shapesetCache.insert(
                typename ShapesetCache::value_type(&shapeset, data)).first;
        }
        const ShapesetData<BasisFunctionType> &shapesetData = shapesetIt->second;

        geometry.getData(Fiber::INTEGRATION_ELEMENTS, shapesetData.quadPoints, geomData);
        for (int pointIndex = 0; pointIndex < shapesetData.quadPoints.n_cols; ++pointIndex)
        {
            for (int functionIndex = 0; 
                 functionIndex < shapesetData.functionCount; ++functionIndex)
                {
                    for (int dim = 0; dim < codomainDim; ++dim)
                    {
                        integral(dim) += 
                            localCoeffs[functionIndex] * 
                            shapesetData.basisData.values(dim, functionIndex, pointIndex) *
                            geomData.integrationElements(pointIndex) *
                            shapesetData.quadWeights[pointIndex];
                    }
                }
        }
    }

    return integral;
}

#define INSTANTIATE_integrateGridFunction(BASIS, RESULT) \
    template \
        arma::Col<RESULT> integrateGridFunction(\
            const GridFunction<BASIS, RESULT>& gridFunction)

FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_integrateGridFunction);

} // end namespace Bempp
