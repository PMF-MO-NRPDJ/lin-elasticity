#pragma once

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
//#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh> // added
//#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <string>
#include <stdexcept>

#include "bctype.hh"
#include "operator.hh"

template <typename GV>
void driver(GV & gv, double E, double nu, double g_vert, double rho, std::string  name)
{
  using namespace Dune::PDELab;

  const int dim = GV::Grid::dimension;
  const int k = 2; // stupanj prostora KE

  // skalarni prostor konačnih elemenata
  using FEM0 = QkLocalFiniteElementMap<GV, double, double, k>;
  using CON  = ConformingDirichletConstraints;
  using VEB0 = ISTL::VectorBackend<>;
  using GFS0 = GridFunctionSpace<GV, FEM0, CON, VEB0>;

  // U vektorskom slučaju dajemo način grupiranja varijabli. Fixed znači da
  // grupiramo komponente koje pripadaju istoj nodalnoj točki.
  using VEB = ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;
  // Vektorski grid function space. Blok varijabli odgovara entitetu! (kod nas vrhu elementa)
  // EntityBlockedOrderingTag = Indicate blocking of the unknowns by grid entity.
  // LexicographicOrderingTag = Indicate lexicographic ordering of the unknowns of non-leaf grid function spaces.
  // InterleavedOrderingTag = Indicate interleaved ordering of the unknowns of non-leaf grid function spaces ...
  using GFS = PowerGridFunctionSpace<GFS0, dim, VEB, Dune::PDELab::EntityBlockedOrderingTag>;

  using CC = typename GFS::template ConstraintsContainer<double>::Type;
  // vektorski rubni uvjeti -- svaka varijabla zadovoljava Dirichletov uvjet na istom dijelu granice.
  using U_BCTypeParam = PowerConstraintsParameters<BCTypeParam<GV>, dim>;

  // lokalni operator
  using LOP = ElasticityLocalOperator<BCTypeParam<GV>, FEM0>;
  using MBE = ISTL::BCRSMatrixBackend<>;
  // konstrukcija vektorskog grid  operatora
  using GO = GridOperator<GFS, GFS, LOP, MBE, double, double, double, CC, CC>;

  // Interpoliramo rubni uvjet
  using BCE0 = BCExtension<GV>;
  // Konstruiraj vektorsku funkciju rubnog uvjeta
  using BCE = PowerGridFunction<BCE0, dim>;
  // Linear solver -- sustav je (skoro) simetričan, možemo koristiti CG_ILU0

  using LS = ISTLBackend_SEQ_CG_ILU0;
   // vektor komponenti
  using U = typename GO::Traits::Domain;
  // linearni solver
  using SLP = StationaryLinearProblemSolver<GO, LS, U>;

  FEM0 fem0(gv);
  GFS0 gfs0(gv,fem0);
  GFS  gfs(gfs0);
  gfs.template child<0>().name("u0");  // Imena nam trebaju za VTK ispis.
  using namespace Dune::Indices;
  //gfs.child(_0).name("u0");
  gfs.child(_1).name("u1");
  gfs.child(_2).name("u2");

  // rubni uvjet za komponentu
  BCTypeParam<GV> bc0(gv);
  U_BCTypeParam bc(bc0);
  // odredi Dirichletovu granicu
  CC cc;
  Dune::PDELab::constraints(bc, gfs, cc);

  // Parametri za lokalni operator
  double mu = E/( 2*(1+nu) );
  double lambda = E*nu/( (1+nu)*(1-2*nu) );
  LOP lop(bc0, mu, lambda, g_vert, rho);
  MBE mbe(std::pow(1 + 2 * k, dim));
  GO go(gfs, cc, gfs, cc, lop, mbe);

  U u(gfs, 0.0);
  BCE0 bce0(gv);   // Dirichletov
  BCE bce(bce0);   // rubni uvjet
  Dune::PDELab::interpolate(bce, gfs, u);

  // ILI ako razne komponente rješenja imaju različite Dirichletove vrijednosti
  // using BCE0 = BCExtension0<GV, double>;
  // using BCE1 = BCExtension0<GV, double>;
  // BCE0 bce0(gv);
  // BCE0 bce1(gv);
  // using BCE = Dune::PDELab::CompositeGridFunction<BCE0, BCE0>;
  // BCE bce(bce0, bce1);
  // Dune::PDELab::interpolate(bce, gfs, u);


  LS ls(5000, true);
  SLP slp(go, ls, u, 1e-8);
  slp.apply();

  if(slp.ls_result().converged)
      std::cout << "Problem solved.\n";
  else
      std::cout << "Solver did not converge.\n";

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, Dune::RefinementIntervals{k});
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, u);  // jednostavniji nain ispisa

  vtkwriter.write(name, Dune::VTK::ascii);

}
