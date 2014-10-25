bempp-simple-vector-spaces
==========================

This repository contains code written in response to issue #168 in the bempp
repository (https://github.com/bempp/bempp/issues/168).

It contains 

* subclasses of Bempp::Space representing spaces of piecewise constant
  and linear vector-valued functions,

* subclasses of Fiber::Shapeset representing bases of constant and linear
  vector-valued functions defined on a reference element and

* a functor class SimpleVectorFunctionValueFunctor that can be used in operators
  acting on vector-valued functions.

The number of vector components of the basis functions is configurable with the
template parameter codomainDim, but at present the templates are explicitly
instantiated only for codomainDim == 3, as at present BEM++ can handle only 3D
geometries. Each basis function has only a single non-zero Cartesian
component. The degrees of freedom are numbered so that the basis function with
index n varies in space as the basis function with index (n / codomainDim) of
the corresponding scalar space, and it is oriented along the (n % codomainDim)th
axis.

Counterparts to these spaces defined on barycentrically refined grid are not
available yet.