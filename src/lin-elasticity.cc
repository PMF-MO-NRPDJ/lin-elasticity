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

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

 // Pročitaj ulaznu datoteku
    Dune::ParameterTree input_data;
    std::string filename (std::string(argv[0])+".input");

    if (argc > 1)
        filename = argv[1];
    try
    {
        Dune::ParameterTreeParser::readINITree (filename, input_data);
    }
    catch (...)
    {
        std::cerr << "The configuration file \"" << filename << "\" "
                  "could not be read. Exiting..." << std::endl;
        std::exit(1);
    }

    int   level   =  input_data.get<int>("level");  // nivo profinjenja
    double E      =  input_data.get<double>("E");   // Youngov modul
    double nu     =  input_data.get<double>("nu");  // Poissonov omjer
    double g_vert =  input_data.get<double>("g_vert");// Površinska sila na presjek
    double rho    =  input_data.get<double>("rho");  // gustoća mase
    std::string name = input_data.get<std::string>("output"); 


    constexpr int dim = 2;  // dimenzija mreže
    using GridType = Dune::YaspGrid<dim>;
    Dune::FieldVector<GridType::ctype,dim> L(2.0);             // Duljina stranice
    L[0] = 20.0;
    std::array<int,dim> s={100,10};          // broj ćelija po stranici
    GridType grid(L, s);  // 20x1
    if(level > 0)
         grid.globalRefine(level);

    auto gv = grid.leafGridView();
    driver(gv, E, nu, g_vert, rho, name);

    return 0;
}
