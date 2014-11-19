%module integrate_grid_function
%{
#define SWIG_FILE_WITH_INIT
#include <numpy/arrayobject.h>
#include "integrate_grid_function.hpp"
%}

%include "bempp.swg"

%init %{
    import_array();
%}

%inline %{
namespace Bempp
{
template <typename BasisFunctionType, typename ResultType>
void _integrateGridFunction(
        const GridFunction<BasisFunctionType, ResultType>& gridFunction,
        arma::Col<ResultType> &result)
{
    result = integrateGridFunction(gridFunction);
}

template <typename BasisFunctionType, typename ResultType>
void _integrateGridFunctionOnSegment(
        const GridFunction<BasisFunctionType, ResultType>& gridFunction,
        const GridSegment& gridSegment,
        arma::Col<ResultType> &result)
{
    result = integrateGridFunctionOnSegment(gridFunction, gridSegment);
}
}
%}

namespace Bempp
{

%apply arma::Col<float>& ARGOUT_COL { arma::Col<float>& result };
%apply arma::Col<double>& ARGOUT_COL { arma::Col<double>& result };
%apply arma::Col<std::complex<float> >& ARGOUT_COL { arma::Col<std::complex<float> >& result };
%apply arma::Col<std::complex<double> >& ARGOUT_COL { arma::Col<std::complex<double> >& result };

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(_integrateGridFunction);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(_integrateGridFunctionOnSegment);

%clear arma::Col<float>& result;
%clear arma::Col<double>& result;
%clear arma::Col<std::complex<float> >& result;
%clear arma::Col<std::complex<double> >& result;
}

%pythoncode %{
    def integrateGridFunction(gridFunction):
        import bempp.lib
        basisFunctionType = gridFunction.basisFunctionType()
        resultType = gridFunction.resultType()
        fullName = ("_integrateGridFunction_" +
                    bempp.lib.checkType(basisFunctionType) + "_" +
                    bempp.lib.checkType(resultType))
        try:
            func = globals()[fullName]
        except KeyError:
            raise TypeError("Function " + fullName + " does not exist.")
        return func(gridFunction)

    def integrateGridFunctionOnSegment(gridFunction, gridSegment):
        import bempp.lib
        basisFunctionType = gridFunction.basisFunctionType()
        resultType = gridFunction.resultType()
        fullName = ("_integrateGridFunctionOnSegment_" +
                    bempp.lib.checkType(basisFunctionType) + "_" +
                    bempp.lib.checkType(resultType))
        try:
            func = globals()[fullName]
        except KeyError:
            raise TypeError("Function " + fullName + " does not exist.")
        return func(gridFunction, gridSegment)
%}


