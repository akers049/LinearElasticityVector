/*
 * NeoHookean_Newton_CompressedStrip.h
 *
 *  Created on: Aug 10, 2017
 *      Author: andrew
 */

#ifndef LINEAR_ELASTIC_H_
#define LINEAR_ELASTIC_H_

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/arpack_solver.h>



#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/slepc_solver.h>


#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Constituitive.h"


#define MAXLINE 1024
#define DIM 2

namespace linear_elastic
{
  using namespace dealii;
  /****************************************************************
                       Class Declarations
  ****************************************************************/


  /****  ElasticProblem  *****
   * This is the primary class used, with all the dealii stuff
   */
  class ElasticProblem
  {
  public:
    ElasticProblem();
    ~ElasticProblem();

    void create_mesh();
    void setup_system ();

    void solve_forward_problem();

    void output_results(const unsigned int cycle) const;

    void read_input_file(char* filename);

    unsigned int get_n_dofs(){return dof_handler.n_dofs();};
    unsigned int get_number_active_cells(){return triangulation.n_active_cells();};


    Vector<double>       present_solution;





  private:

    void setup_system_constraints();
    void solve_adjoint_problem();

    void apply_boundaries_and_constraints();

    void assemble_system_matrix();
    void assemble_system_rhs();

    void solve();

    void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
        int const maxSize, int* const endOfFileFlag);



    Triangulation<DIM,DIM>   triangulation;
    DoFHandler<DIM>      dof_handler;

    FESystem<DIM>        fe;

    ConstraintMatrix     constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       system_rhs;

    char output_directory[MAXLINE];

    std::vector<unsigned int>  grid_dimensions;
    std::vector<double> domain_dimensions;

    LinearElastic *LE = NULL;



    double E = 1.0;
    double nu = 1.0;
    unsigned int corner_dof = 0;

    unsigned int qx = 2;
    unsigned int qy = 2;

  };
}

#endif /* LINEAR_ELASTIC_H_ */
