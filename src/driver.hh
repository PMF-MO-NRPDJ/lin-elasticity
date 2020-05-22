#ifndef DRIVER_HH_IS_INCLUDED
#define DRIVER_HH_IS_INCLUDED

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
  const int dim = GV::Grid::dimension;
  const int k = 2; // stupanj prostora KE
  // skalarni prostor konačnih elemenata
  using FEM0 = Dune::PDELab::QkLocalFiniteElementMap<GV, double, double, k>;
  //using CON = Dune::PDELab::NoConstraints;
  using CON = Dune::PDELab::ConformingDirichletConstraints;
  using VEB0 = Dune::PDELab::ISTL::VectorBackend<>;
  using GFS0 = Dune::PDELab::GridFunctionSpace<GV, FEM0, CON, VEB0>;
  // U vektorskom slučaju dajemo način grupiranja varijabli. Fixed znači da
  // grupiramo komponente koje pripadaju istoj nodalnoj točki.
  using VEB = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;
  // Vektorski grid function space. Blok varijabli odgovara entitetu! (kod nas vrhu elementa)
  // EntityBlockedOrderingTag = Indicate blocking of the unknowns by grid entity.
  // LexicographicOrderingTag = Indicate lexicographic ordering of the unknowns of non-leaf grid function spaces.
  // InterleavedOrderingTag = Indicate interleaved ordering of the unknowns of non-leaf grid function spaces ...
  using GFS = Dune::PDELab::PowerGridFunctionSpace<GFS0, dim, VEB, Dune::PDELab::EntityBlockedOrderingTag>;
  using CC = typename GFS::template ConstraintsContainer<double>::Type;
  // vektorski rubni uvjeti -- svaka varijabla zadovoljava Dirichletov uvjet na istom dijelu granice.
  using U_BCTypeParam = Dune::PDELab::PowerConstraintsParameters<BCTypeParam<GV>, dim>;
  // lokalni operator
  using LOP = ElasticityLocalOperator<BCTypeParam<GV>, FEM0>;
  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  // konstrukcija vektorskog grid  operatora
  using GO = Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, double, double, double, CC, CC>;
  // Interpoliramo rubni uvjet
  using BCE0 = BCExtension<GV, double>;
  // Konstruiraj vektorsku funkciju rubnog uvjeta
  using BCE = Dune::PDELab::PowerGridFunction<BCE0, dim>;
  // Linear solver -- sustav je (skoro) simetričan, možemo koristiti CG_ILU0
  using LS = Dune::PDELab::ISTLBackend_SEQ_CG_ILU0;
   // vektor komponenti
  using U = typename GO::Traits::Domain;
  // linearni solver
  using SLP = Dune::PDELab::StationaryLinearProblemSolver<GO, LS, U>;
  // Uzmimo tipove za potprostore: 0 -> potprostor prve komponente, 1 -> druge itd.
  using U0SUB = Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<0> >;
  using U1SUB = Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<1> >;
  // Napravi mrežnu funciju od svake komponenete rješenja
  using U0_DGF =  Dune::PDELab::DiscreteGridFunction<U0SUB, U> ;
  using U1_DGF =  Dune::PDELab::DiscreteGridFunction<U1SUB, U> ;

  FEM0 fem0(gv);
  GFS0 gfs0(gv,fem0);
  GFS  gfs(gfs0);
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

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, Dune::RefinementIntervals{1});
  U0SUB u0sub(gfs); // prostor za prvu komponentu
  U1SUB u1sub(gfs); // prostor za drugu komponentu
  U0_DGF u0_dgf(u0sub, u);
  U1_DGF u1_dgf(u1sub, u);

  // Ispiši mrežne funkcije
  vtkwriter.addVertexData( std::make_unique<Dune::PDELab::VTKGridFunctionAdapter<U0_DGF>>(u0_dgf, "u_x"));
  vtkwriter.addVertexData( std::make_unique<Dune::PDELab::VTKGridFunctionAdapter<U1_DGF>>(u1_dgf, "u_y"));
  vtkwriter.write(name, Dune::VTK::ascii);

}
#endif   
