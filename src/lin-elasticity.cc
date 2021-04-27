#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include "driver.hh"
/* Problem linearizirane eleastičnosti.
 * Domena je greda (0,20)x(0,1)x(0,2) u polju sile teže s
 * učvršćena u X=0, s opterećenjem na gornjoj stranici.
 * Računamo pomak.
 */
int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    // Pročitaj ulaznu datoteku
    Dune::ParameterTree input_data;
    std::string filename (std::string(argv[0])+".input");
    Dune::ParameterTreeParser::readINITree (filename, input_data);

    // Pročitaj parametre
    int         level  = input_data.get<int>        ("level");   // nivo profinjenja
    double      E      = input_data.get<double>     ("E");       // Youngov modul
    double      nu     = input_data.get<double>     ("nu");      // Poissonov omjer
    double      g_vert = input_data.get<double>     ("g_vert");  // Površinska sila na presjek
    double      rho    = input_data.get<double>     ("rho");     // gustoća mase
    std::string name   = input_data.get<std::string>("output");  // ime izlazne datoteke

    // Konstruiraj mrežu
    constexpr int dim = 3;  // dimenzija mreže
    using GridType = Dune::YaspGrid<dim>;
    Dune::FieldVector<GridType::ctype,dim> L;             // Duljina stranice
    L[0] = 20.0;
    L[1] = 1.0;
    L[2] = 2.0;
    std::array<int,dim> s={100,5,10};          // broj ćelija po stranici
    GridType grid(L, s);
    if(level > 0)
         grid.globalRefine(level);

    // Zovi driver.
    auto gv = grid.leafGridView();
    driver(gv, E, nu, g_vert, rho, name);

    return 0;
}
