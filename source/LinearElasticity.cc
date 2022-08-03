/* ---------------------------------------------------------------------
 *
 *
 *
 *
 *
 * ---------------------------------------------------------------------
 *
 *
 * Author: Andrew Akerson
 */

#ifndef LINEAR_ELASTIC_CC_
#define LINEAR_ELASTIC_CC_
#include "LinearElasticity.h"

#include <fstream>
#include <iostream>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <bits/stdc++.h>


#define DIM 2



namespace linear_elastic
{
  using namespace dealii;

  ElasticProblem::ElasticProblem ()
    :
    dof_handler (triangulation),
    fe (FESystem<DIM>(FE_Q<DIM>(1), DIM), 1)
  {}




  ElasticProblem::~ElasticProblem ()
  {
    dof_handler.clear ();

    delete LE;

  }


  void ElasticProblem::create_mesh()
  {

    // creates our strip.
    Point<DIM> corner1, corner2;
    corner1(0) =  0.0;
    corner1(1) =  0.0;
    corner2(0) =  domain_dimensions[0];
    corner2(1) =  domain_dimensions[1];

    GridGenerator::subdivided_hyper_rectangle (triangulation, grid_dimensions, corner1, corner2, true);

  }

  void ElasticProblem::setup_system ()
  {
    // Sets up system. Makes the constraint matrix, and reinitializes the
    // vectors and matricies used throughout to be the proper size.


    dof_handler.distribute_dofs (fe);

    present_solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    setup_system_constraints();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);



    sparsity_pattern.copy_from (dsp);


    system_matrix.reinit (sparsity_pattern);

    LE = new LinearElastic(E, nu);



  }


  void ElasticProblem::setup_system_constraints ()
  {

    constraints.clear ();

    constraints.close ();

    // now do hanging nodes. Because some of the constraints might refer to the same dof
    // for both the symmetry constraint and the hanging node constraint, we will make them
    // separate, then merge them, giving precedence to the hanging node constraints;
    ConstraintMatrix hanging_node_constraints;
    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
    hanging_node_constraints.close();

    constraints.merge(hanging_node_constraints, ConstraintMatrix::MergeConflictBehavior::right_object_wins);
  }




  void ElasticProblem::solve_forward_problem()
  {

    assemble_system_matrix();
    assemble_system_rhs();
    apply_boundaries_and_constraints();
    solve();
    output_results(42069);
  }

  void ElasticProblem::assemble_system_matrix()
  {
    // Assembling the system matrix. I chose to make the rhs and system matrix assemblies separate,
    // because we only do one at a time anyways in the newton method.

    system_matrix = 0.0;

    QGauss<1> quad_x(qx);
    QGauss<1> quad_y(qy);


    QAnisotropic<DIM> quadrature_formula(quad_x, quad_y);


    FEValues<DIM> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;

    unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const FEValuesExtractors::Vector u(0);

    typename DoFHandler<DIM>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {

      cell_matrix = 0.0;

      fe_values.reinit (cell);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        for (unsigned int n = 0; n < dofs_per_cell; ++n)
          for (unsigned int m = 0; m < dofs_per_cell; ++m)
          {

            for(unsigned int i = 0; i < DIM; i ++)
              for(unsigned int j = 0; j < DIM; j ++)
                  for(unsigned int k = 0; k < DIM; k ++)
                    for (unsigned int l = 0; l<DIM; l ++)
                    {

                      cell_matrix(n,m) += LE->C[i][j][k][l]*fe_values[u].symmetric_gradient(n, q_point)[i][j]
                                                           *fe_values[u].symmetric_gradient(m, q_point)[k][l]*fe_values.JxW(q_point);
                    }
          }

      }

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int n=0; n<dofs_per_cell; ++n)
        for (unsigned int m=0; m<dofs_per_cell; ++m)
        {
          system_matrix.add (local_dof_indices[n],
                             local_dof_indices[m],
                             cell_matrix(n,m));
        }
    }



  }

  void ElasticProblem::assemble_system_rhs()
  {
    Point<DIM> corner_point(domain_dimensions[0], 0.0);
    Point<DIM> dir(0.0, -0.01);

    VectorTools::create_point_source_vector(dof_handler, corner_point, dir, system_rhs);
  }

  void ElasticProblem::apply_boundaries_and_constraints()
  {
    constraints.condense (system_matrix);
    constraints.condense (system_rhs);


    std::map<types::global_dof_index,double> boundary_values;

    std::vector<bool> side2_components = {true, true};
    ComponentMask side2_mask(side2_components);

    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              ZeroFunction<DIM, double>(DIM),
                                              boundary_values,
                                              side2_mask);

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        present_solution,
                                        system_rhs);
  }


  void ElasticProblem::solve ()
  {


    // direct solver for the system
    if(dof_handler.n_dofs() < 10000)
    {
      SparseDirectUMFPACK  A_direct;
      A_direct.initialize(system_matrix);
      A_direct.vmult (present_solution, system_rhs);
    }
    else
    {
      SolverControl solver_control(dof_handler.n_dofs(), 1e-11);
      SolverCG<> solver(solver_control);
      solver.solve(system_matrix, present_solution, system_rhs, PreconditionIdentity());
    }

    constraints.distribute (present_solution);


  }



  void ElasticProblem::output_results (const unsigned int cycle) const
  {

    std::vector<std::string> solution_names;
    switch (DIM)
      {
      case 1:
        solution_names.push_back ("displacement");
        break;
      case 2:
        solution_names.push_back ("x1_displacement");
        solution_names.push_back ("x2_displacement");
        break;
      case 3:
        solution_names.push_back ("x1_displacement");
        solution_names.push_back ("x2_displacement");
        solution_names.push_back ("x3_displacement");
        break;
      default:
        Assert (false, ExcNotImplemented());
        break;
      }


    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(DIM,
                     DataComponentInterpretation::component_is_part_of_vector);

    std::string filename0(output_directory);
    filename0 += "/lagrangian_solution";

    // see if the directory exists...
    struct stat st;
    if (stat(filename0.c_str(), &st) == -1)
      mkdir(filename0.c_str(), 0700);

    filename0 += "/lagrangian_solution-";
    filename0 += std::to_string(cycle);
    filename0 += ".vtk";
    std::ofstream output_lagrangian_solution (filename0.c_str());

    DataOut<DIM> data_out_lagrangian;


    data_out_lagrangian.add_data_vector(dof_handler,
                                       present_solution,
                                       solution_names,
                                       interpretation);

    data_out_lagrangian.build_patches ();
    data_out_lagrangian.write_vtk (output_lagrangian_solution);

  }



  void ElasticProblem::read_input_file(char* filename)
  {
    FILE* fid;
    int endOfFileFlag;
    char nextLine[MAXLINE];

    int valuesWritten;
    bool fileReadErrorFlag = false;

    grid_dimensions.resize(DIM);
    domain_dimensions.resize(DIM);

    fid = std::fopen(filename, "r");
    if (fid == NULL)
    {
      std::cout << "Unable to open file \"" << filename  << "\"" <<  std::endl;
      fileReadErrorFlag = true;
    }
    else
    {

      // Read in the output name
      char directory_name[MAXLINE];
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%s", directory_name);
      if (valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      sprintf(output_directory, "output/");
      strcat(output_directory, directory_name);

      // Read in the grid dimensions
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u %u", &grid_dimensions[0], &grid_dimensions[1]);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // Read in the domain dimensions
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg", &domain_dimensions[0], &domain_dimensions[1]);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // Read in E and nu
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg", &E, &nu);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the number of guass points in the x and y direction
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u  %u", &qx, &qy);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      fileClose:
      {
        fclose(fid);
      }
    }

    if (fileReadErrorFlag)
    {
      // default parameter values
      std::cout << "Error reading input file, Exiting.\n" << std::endl;
      exit(1);
    }
    else
      std::cout << "Input file successfully read" << std::endl;

    // make the output directory
    struct stat st;
    if (stat("./output", &st) == -1)
       mkdir("./output", 0700);

    if (stat(output_directory, &st) == -1)
      mkdir(output_directory, 0700);

  }

  void ElasticProblem::getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                          int const maxSize, int* const endOfFileFlag)
  {
    *endOfFileFlag = 0;
    do
    {
      if(fgets(nextLinePtr, maxSize, filePtr) == NULL)
      {
        *endOfFileFlag = 1;
        break;
      }
      while ((nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t') ||
             (nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r' ))
      {
        nextLinePtr = (nextLinePtr + 1);
      }
    }
    while ((strncmp("#", nextLinePtr, 1) == 0) || (strlen(nextLinePtr) == 0));
  }

}

#endif // LINEAR_ELASTIC_CC_
