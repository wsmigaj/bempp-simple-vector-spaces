%module simple_vector_spaces
%{
#define SWIG_FILE_WITH_INIT
#include <numpy/arrayobject.h>
#include "piecewise_constant_vector_space.hpp"
#include "piecewise_linear_continuous_vector_space.hpp"
#include "piecewise_linear_discontinuous_vector_space.hpp"
%}

%include "bempp.swg"

%inline %{
namespace Bempp
{

    template <typename BasisFunctionType>
        boost::shared_ptr< Space< BasisFunctionType > >
        piecewiseConstantVectorSpace(
            const boost::shared_ptr<const Grid>& grid,
            const GridSegment* segment = 0)
    {
        typedef PiecewiseConstantVectorSpace<BasisFunctionType, 3> Type;
        if (segment)
            return boost::shared_ptr<Type>(new Type(grid, *segment));
        else
            return boost::shared_ptr<Type>(new Type(grid));
    }

    template <typename BasisFunctionType>
        boost::shared_ptr< Space< BasisFunctionType > >
        piecewiseLinearContinuousVectorSpace(
            const boost::shared_ptr<const Grid>& grid,
            const GridSegment* segment = 0,
            bool strictlyOnSegment = false)
    {
        typedef PiecewiseLinearContinuousVectorSpace<BasisFunctionType, 3> Type;
        if (segment)
            return boost::shared_ptr<Type>(new Type(grid, *segment, strictlyOnSegment));
        else
            return boost::shared_ptr<Type>(new Type(grid));
    }

    template <typename BasisFunctionType>
        boost::shared_ptr< Space< BasisFunctionType > >
        piecewiseLinearDiscontinuousVectorSpace(
            const boost::shared_ptr<const Grid>& grid,
            const GridSegment* segment = 0,
            bool strictlyOnSegment = false)
    {
        typedef PiecewiseLinearDiscontinuousVectorSpace<BasisFunctionType, 3> Type;
        if (segment)
            return boost::shared_ptr<Type>(new Type(grid, *segment, strictlyOnSegment));
        else
            return boost::shared_ptr<Type>(new Type(grid));
    }
}
%}

%feature("compactdefaultargs") Bempp::piecewiseConstantVectorSpace;
%feature("compactdefaultargs") Bempp::piecewiseLinearContinuousVectorSpace;
%feature("compactdefaultargs") Bempp::piecewiseLinearDiscontinuousVectorSpace;
%template(createPiecewiseConstantVectorSpace) Bempp::piecewiseConstantVectorSpace<double>;
%template(createPiecewiseLinearContinuousVectorSpace) Bempp::piecewiseLinearContinuousVectorSpace<double>;
%template(createPiecewiseLinearDiscontinuousVectorSpace) Bempp::piecewiseLinearDiscontinuousVectorSpace<double>;

