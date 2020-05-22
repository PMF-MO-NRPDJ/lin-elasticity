add_executable("lin-elasticity" lin-elasticity.cc driver.hh operator.hh bctype.hh lin-elasticity.input)
target_link_dune_default_libraries("lin-elasticity")

dune_symlink_to_source_files(FILES lin-elasticity.input)
