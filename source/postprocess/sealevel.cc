/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/sealevel.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/simulator.h>
#include <aspect/global.h>
#include <aspect/postprocess/geoid.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <cmath>
#include <limits>

// TODO:
// Obtain geoid anomaly of surface load
// Compute new surface gravity per location, dependent on surface topography and Earth model

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    void
    SeaLevel<dim>::initialize()
    {
      topography_lookup = std::make_unique<Utilities::StructuredDataLookup<dim-1>>(1,1.0);
      topography_lookup->load_file(data_directory_topography+data_file_name_topography,this->get_mpi_communicator());

      iceheight_lookup = std::make_unique<Utilities::StructuredDataLookup<dim-1>>(1,1.0);
      iceheight_lookup->load_file(data_directory_iceheight+data_file_name_iceheight,this->get_mpi_communicator());
    }


    template <int dim>
    double
    SeaLevel<dim>::compute_nonuniform_sealevel_change(const Point<dim> &position) const
    {
      const double geoid_displacement = Geoid<dim>::evaluate(position); // (check sign of geoid_displacement)  
      const double elevation = this->get_geometry_model().height_above_reference_surface(position);
      
      Point<dim-1> internal_position;
      const std::array<double,dim> spherical_position = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
      for (unsigned int d = 1; d < dim; ++d)
        internal_position[d-1] = spherical_position[d];

      topography = topography_lookup->get_data(internal_position,0);
      if (topography > 0.0)
        const double ocean_mask = 1;
      else
        const double ocean_mask = 0;
      
      nonuniform_sealevel = (geoid_displacement-elevation)*ocean_mask;

      return nonuniform_sealevel_change;
    }


    template <int dim>
    double
    SeaLevel<dim>::compute_sealevel_offset(const double &outer_radius) const
    {
      const unsigned int quadrature_degree = this->introspection().polynomial_degree.temperature;
      const QGauss<2> quadrature_formula_face(quadrature_degree);

      FEFaceValues<3> fe_face_values (this->get_mapping(),
                                      this->get_fe(),
                                      quadrature_formula_face,
                                      update_values |
                                      update_quadrature_points |
                                      update_JxW_values);


      // Vectors to store the location, infinitesimal area, and topography/geoid displacement/ocean mask/ice height associated with each quadrature point of each surface cell.
      std::vector<std::pair<Point<3>,std::pair<double,double>>> topo_stored_values;
      std::vector<std::pair<Point<3>,std::pair<double,double>>> geoid_stored_values;
      std::vector<std::pair<Point<3>,std::pair<double,double>>> oceanmask_stored_values;
      std::vector<std::pair<Point<3>,std::pair<double,double>>> iceheight_stored_values;

      // Loop over all of the boundary cells and if one is at the
      // surface, evaluate the topography/geoid displacement/ocean mask/ice height there.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          {
            bool at_upper_surface = false;
            {
              for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
                {
                  if (cell->at_boundary(f) && cell->face(f)->boundary_id() == this->get_geometry_model().translate_symbolic_boundary_name_to_id("top"))
                    {
                      at_upper_surface = true;
                      break;
                    }
                  else
                    continue;
                }

              // If the cell is not at the top boundary, jump to the next cell.
              if (at_upper_surface == false)
                continue;
            }

            // Focus on the boundary cell's upper face if on the top boundary.
            fe_face_values.reinit(cell);

            // If the cell is at the top boundary, add its contributions to the topography/geoid displacement/ocean mask/ice height storage vectors.
            if (at_upper_surface)
              {
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    const Point<3> current_position = fe_face_values.quadrature_point(q);
                    const double topography = this->get_geometry_model().height_above_reference_surface(current_position);
                    topo_stored_values.emplace_back (current_position, std::make_pair(fe_face_values.JxW(q), topography));
                    const double geoid = Geoid<dim>::evaluate(current_position);
                    geoid_stored_values.emplace_back (fe_face_values.quadrature_point(q), std::make_pair(fe_face_values.JxW(q), geoid);

                    Point<dim-1> internal_position;
                    const std::array<double,dim> spherical_position = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
                    for (unsigned int d = 1; d < dim; ++d)
                      internal_position[d-1] = spherical_position[d];

                    topography = topography_lookup->get_data(internal_position,0);
                    if (topography > 0.0)
                      const double ocean_mask = 1;
                    else
                      const double ocean_mask = 0;
                    oceanmask_stored_values.emplace_back (fe_face_values.quadrature_point(q), std::make_pair(fe_face_values.JxW(q), ocean_mask);

                    ice_height = iceheight_lookup->get_data(internal_position,0);
                    iceheight_stored_values.emplace_back (fe_face_values.quadrature_point(q), std::make_pair(fe_face_values.JxW(q), ice_height);
                    
                    // Compute required integrals for the sea level offset.
                    const double integral_oceanmask += ocean_mask*fe_face_values.JxW(q);
                    const double integral_iceheight += iceheight*density_ice*(1-ocean_mask)*fe_face_values.JxW(q);
                    const double integral_topo_geoid += (geoid-topography)*ocean_mask*fe_face_values.JxW(q);
                  }
              }
          }
      
      const double sealevel_offset = -1./integral_oceanmask*(1./density_water*integral_iceheight+integral_topo_geoid);

      std::vector<std::vector<double>> topo_spherical_function;
      std::vector<std::vector<double>> geoid_spherical_function;
      std::vector<std::vector<double>> oceanmask_spherical_function;
      std::vector<std::vector<double>> iceheight_spherical_function;

      for (unsigned int i=0; i<topo_stored_values.size(); ++i)
        {
          const std::array<double,3> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(topo_stored_values[i].first);

          // Calculate spherical infinitesimal sin(theta)*d_theta*d_phi by infinitesimal_area/radius^2
          const double infinitesimal = topo_stored_values[i].second.first/(outer_radius*outer_radius);

          // Theta, phi, spherical infinitesimal, and topography
          topo_spherical_function.emplace_back(std::vector<double> {scoord[2],
                                                                    scoord[1],
                                                                    infinitesimal,
                                                                    topo_stored_values[i].second.second
                                                                   });
          
          // Theta, phi, spherical infinitesimal, and geoid displacement  
          geoid_spherical_function.emplace_back(std::vector<double> {scoord[2],
                                                                     scoord[1],
                                                                     infinitesimal,
                                                                     geoid_stored_values[i].second.second
                                                                    });

          // Theta, phi, spherical infinitesimal, and ocean mask  
          oceanmask_spherical_function.emplace_back(std::vector<double> {scoord[2],
                                                                         scoord[1],
                                                                         infinitesimal,
                                                                         oceanmask_stored_values[i].second.second
                                                                        });

          // Theta, phi, spherical infinitesimal, and ice height
          iceheight_spherical_function.emplace_back(std::vector<double> {scoord[2],
                                                                         scoord[1],
                                                                         infinitesimal,
                                                                         iceheight_stored_values[i].second.second
                                                                        });
        }

      return sealevel_offset
    }


    template <int dim>
    double
    SeaLevel<dim>::compute_total_surface_pressure(const Point<dim> &position) const
    {
      const double nonuniform_sealevel_change = compute_nonuniform_sealevel_change(position);
      
      Point<dim-1> internal_position;
      const std::array<double,dim> spherical_position = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
      for (unsigned int d = 1; d < dim; ++d)
        internal_position[d-1] = spherical_position[d];
      ice_height = iceheight_lookup->get_data(internal_position,0);

      Point<dim> surface_point;
      surface_point[0] = outer_radius;
      const double surface_gravity = this->get_gravity_model().gravity_vector(surface_point).norm();

      const double total_surface_pressure = -surface_gravity*((nonuniform_sealevel_change+sealevel_offset)*density_water+ice_height*density_ice);

      return total_surface_pressure
    }


    template <int dim>
    std::pair<std::string,std::string>
    SeaLevel<dim>::execute (TableHandler &statistics)
    {
      // Disallow use of the plugin for any other than a 3D spherical shell geometry.
      AssertThrow (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model())
                   &&
                   dim == 3,
                   ExcMessage("The sea level postprocessor is currently only implemented for the 3D spherical shell geometry model."));

      const GeometryModel::SphericalShell<dim> &geometry_model = Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model());
      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");

      // Get the value of the outer radius.
      const double outer_radius = geometry_model.outer_radius();
      
      // Get the sea level offset (constant for every location).
      sealevel_offset = compute_sealevel_offset(outer_radius);

      // Writing output data: non-uniform sea level change (add other output data later too).
      // Get a quadrature rule that exists only on the corners.
      QTrapezoid<dim-1> face_corners;
      FEFaceValues<dim> face_vals (this->get_mapping(), this->get_fe(), face_corners, update_quadrature_points);

      // Have a stream into which we write the data. The text stream is then later send to processor 0.
      std::ostringstream output_stats;
      std::ostringstream output_file;

      // Choose stupidly large values for initialization.
      double local_max_height = -std::numeric_limits<double>::max();
      double local_min_height = std::numeric_limits<double>::max();

      // Loop over all of the surface cells and save the non-uniform sea level change to stored_value.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
            if (cell->face(face_no)->at_boundary())
              {
                if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                  continue;

                face_vals.reinit( cell, face_no);

                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                    const Point<dim> vertex = face_vals.quadrature_point(corner);
                    const double nonuniform_sealevel_change = compute_nonuniform_sealevel_change(vertex);
                    if (write_to_file)
                      output_file << vertex << ' '<< elevation << std::endl;
                    if ( nonuniform_sealevel_change > local_max_height)
                      local_max_height = nonuniform_sealevel_change;
                    if ( nonuniform_sealevel_change < local_min_height)
                      local_min_height = nonuniform_sealevel_change;
                  }
              }

      // Calculate min/max non uniform sea level change across all processes.
      const double max_nonuniform_sealevel_change = Utilities::MPI::max(local_max_height, this->get_mpi_communicator());
      const double min_nonuniform_sealevel_change = Utilities::MPI::min(local_min_height, this->get_mpi_communicator());

      // Write results to statistics file.
      statistics.add_value ("Minimum non-uniform sea level change (m)",
                            min_nonuniform_sealevel_change);
      statistics.add_value ("Maximum non-uniform sea level change (m)",
                            max_nonuniform_sealevel_change);
      const char *columns[] = { "Minimum non-uniform sea level change (m)",
                                "Maximum non-uniform sea level change (m)"
                              };
      for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
        {
          statistics.set_precision (columns[i], 8);
          statistics.set_scientific (columns[i], true);
        }

      output_stats.precision(4);
      output_stats << min_nonuniform_sealevel_change << " m, "
                   << max_nonuniform_sealevel_change << " m";

      // Write the solution to file.

      // If this is the first time we get here, set the last output time
      // to the current time - output_interval. This makes sure we
      // always produce data during the first time step.
      if (std::isnan(last_output_time))
        {
          last_output_time = this->get_time() - output_interval;
        }

      // Just return stats if text output is not required at all or not needed at this time.
      if (!write_to_file || ((this->get_time() < last_output_time + output_interval)
                             && (this->get_timestep_number() != 0)))
        return std::pair<std::string, std::string> ("Non-uniform sea level change min/max:",
                                                    output_stats.str());

      std::string filename = this->get_output_directory() +
                             "nonuniform_sealevel_change." +
                             Utilities::int_to_string(this->get_timestep_number(), 5);
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        filename.append("." + Utilities::int_to_string (this->get_nonlinear_iteration(), 4));

      const unsigned int max_data_length = Utilities::MPI::max (output_file.str().size()+1,
                                                                this->get_mpi_communicator());

      const unsigned int mpi_tag = 777;

      // On processor 0, collect all of the data the individual processors sent
      // and concatenate them into one file.
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::ofstream file (filename.c_str());

          file << "# "
               << "x y z"
               << " nonuniform_sealevel_change" << std::endl;

          // First write out the data we have created locally.
          file << output_file.str();

          std::string tmp;
          tmp.resize (max_data_length, '\0');

          // Then loop through all of the other processors and collect
          // data, then write it to the file.
          for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
              MPI_Status status;
              // Get the data. Note that MPI says that an MPI_Recv may receive
              // less data than the length specified here. Since we have already
              // determined the maximal message length, we use this feature here
              // rather than trying to find out the exact message length with
              // a call to MPI_Probe.
              const int ierr = MPI_Recv (&tmp[0], max_data_length, MPI_CHAR, p, mpi_tag,
                                         this->get_mpi_communicator(), &status);
              AssertThrowMPI(ierr);

              // Output the string. Note that 'tmp' has length max_data_length,
              // but we only wrote a certain piece of it in the MPI_Recv, ended
              // by a \0 character. Write only this part by outputting it as a
              // C string object, rather than as a std::string.
              file << tmp.c_str();
            }
        }
      else
        // On other processors, send the data to processor zero. include the \0
        // character at the end of the string.
        {
          output_file << "\0";
          const int ierr = MPI_Send (&output_file.str()[0], output_file.str().size()+1, MPI_CHAR, 0, mpi_tag,
                                     this->get_mpi_communicator());
          AssertThrowMPI(ierr);
        }

      // if output_interval is positive, then update the last supposed output time
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // This is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((this->get_time()-last_output_time)/output_interval*magic) * output_interval/magic;
        }

      return std::pair<std::string, std::string> ("Non-uniform sea level change min/max:",
                                                  output_stats.str());
    }


    template <int dim>
    std::list<std::string>
    SeaLevel<dim>::required_other_postprocessors() const
    {
      std::list<std::string> deps;
      deps.emplace_back("geoid");
      return deps;
    }


    template <int dim>
    double
    Sealevel<dim>::evaluate (const Point<dim> &/*position*/) const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }


    template <int dim>
    void
    SeaLevel<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Sea level");
        {
          prm.declare_entry ("Ice density", "931",
                              Patterns::Anything (),
                             "The density of ice [kg/m3]");
          prm.declare_entry ("Water density", "1000",
                              Patterns::Anything (),
                             "The density of water [kg/m3]");

          prm.declare_entry ("Data directory topography",
                             "$ASPECT_SOURCE_DIR/data/geometry-model/initial-topography-model/ascii-data/test/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the topography "
                             "ascii data. This path may either be absolute (if starting with a "
                             "`/') or relative to the current directory. The path may also "
                             "include the special text `$ASPECT_SOURCE_DIR' which will be "
                             "interpreted as the path in which the ASPECT source files were "
                             "located when ASPECT was compiled. This interpretation allows, "
                             "for example, to reference files located in the `data/' subdirectory "
                             "of ASPECT.");
          prm.declare_entry ("Data file name topography", "shell_3d_outer.0.txt",
                             Patterns::Anything (),
                             "The file name of the topography ascii data. For the ascii data, "
                             "provide file in the same format as described in "
                             "'ascii data' initial composition plugin." );

          prm.declare_entry ("Data directory ice height",
                             "$ASPECT_SOURCE_DIR/data/boundary_traction/ascii-data/test/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the ice height [m] "
                             "ascii data. This path may either be absolute (if starting with a "
                             "`/') or relative to the current directory. The path may also "
                             "include the special text `$ASPECT_SOURCE_DIR' which will be "
                             "interpreted as the path in which the ASPECT source files were "
                             "located when ASPECT was compiled. This interpretation allows, "
                             "for example, to reference files located in the `data/' subdirectory "
                             "of ASPECT.");
          prm.declare_entry ("Data file name ice height", "shell_3d_outer.0.txt",
                             Patterns::Anything (),
                             "The file name of the ice height ascii data. For the ascii data, "
                             "provide file in the same format as described in "
                             "'ascii data' initial composition plugin." );

          prm.declare_entry ("Output to file", "false",
                             Patterns::List(Patterns::Bool()),
                             "Whether or not to write sea level to a text file named named "
                             "'sealevel.NNNNN' in the output directory");
          prm.declare_entry ("Time between text output", "0.",
                             Patterns::Double (0.),
                             "The time interval between each generation of "
                             "text output files. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    SeaLevel<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Sea level");
        {
          density_ice = prm.get ("Ice density");
          density_water = prm.get ("Water density");
          data_directory_topography = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory topography"));
          data_file_name_topography = prm.get ("Data file name topography");
          data_directory_iceheight = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory ide load"));
          data_file_name_iceheight = prm.get ("Data file name ice height");
          write_to_file = prm.get_bool ("Output to file");
          output_interval = prm.get_double ("Time between text output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(SeaLevel,
                                  "sea level",
                                  "A postprocessor intended for use with a deforming top surface. After every step "
                                  "it computes the sea level based on the topography, ocean basin, ice melt, "
                                  "perturbed gravitational potential of the Earth model and gravitational potential "
                                  "of the ice load, relative to a reference datum (initial "
                                  "radius for a spherical shell geometry model). "
                                  "If 'SeaLevel.Output to file' is set to true, also outputs sea level "
                                  "into text files named `sealevel.NNNNN' in the output directory, "
                                  "where NNNNN is the number of the time step.\n"
                                  "The file format then consists of lines with Euclidean coordinates "
                                  "followed by the corresponding sea level value."
                                  "Sea level is printed/written in meters.")
  }
}

      // // Get the non-uniform sea level change.
      // const unsigned int quadrature_degree = this->introspection().polynomial_degree.temperature;
      // const QGauss<2> quadrature_formula_face(quadrature_degree);

      // FEFaceValues<3> fe_face_values (this->get_mapping(),
      //                                 this->get_fe(),
      //                                 quadrature_formula_face,
      //                                 update_values |
      //                                 update_quadrature_points |
      //                                 update_JxW_values);

      // // Vectors to store the location, infinitesimal area, and non-uniform sea level change associated with each quadrature point of each surface cell.
      // std::vector<std::pair<Point<3>,std::pair<double,double>>> nonuniform_sealevel_change_stored_values;

      // // Loop over all of the boundary cells and if one is at the
      // // surface, evaluate the non-uniform sea level change there.
      // for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      //   if (cell->is_locally_owned() && cell->at_boundary())
      //     {
      //       bool at_upper_surface = false;
      //       {
      //         for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
      //           {
      //             if (cell->at_boundary(f) && cell->face(f)->boundary_id() == this->get_geometry_model().translate_symbolic_boundary_name_to_id("top"))
      //               {
      //                 at_upper_surface = true;
      //                 break;
      //               }
      //             else
      //               continue;
      //           }

      //         // If the cell is not at the top boundary, jump to the next cell.
      //         if (at_upper_surface == false)
      //           continue;
      //       }

      //       // Focus on the boundary cell's upper face if on the top boundary.
      //       fe_face_values.reinit(cell);

      //       // If the cell is at the top boundary, add its contributions to the non-uniform sea level change storage vectors.
      //       if (at_upper_surface)
      //         {
      //           for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
      //             {
      //               const Point<3> current_position = fe_face_values.quadrature_point(q);
      //               const double nonuniform_sealevel_change = compute_nonuniform_sealevel_change(current_position);
      //               nonuniform_sealevel_change_stored_values.emplace_back (current_position, std::make_pair(fe_face_values.JxW(q), nonuniform_sealevel_change));
      //             }
      //         }
      //     }
      
      // std::vector<std::vector<double>> nonuniform_sealevel_change_spherical_function;

      // for (unsigned int i=0; i<nonuniform_sealevel_change_stored_values.size(); ++i)
      //   {
      //     const std::array<double,3> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(nonuniform_sealevel_change_stored_values[i].first);

      //     // Calculate spherical infinitesimal sin(theta)*d_theta*d_phi by infinitesimal_area/radius^2
      //     const double infinitesimal = nonuniform_sealevel_change_stored_values[i].second.first/(outer_radius*outer_radius);

      //     // Theta, phi, spherical infinitesimal, and non-uniform sea level change
      //     nonuniform_sealevel_change_spherical_function.emplace_back(std::vector<double> {scoord[2],
      //                                                                                     scoord[1],
      //                                                                                     infinitesimal,
      //                                                                                     nonuniform_sealevel_change_stored_values[i].second.second
      //                                                                                    });
      //   }