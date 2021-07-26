function [F,V] = execute_triangle(command_line_args, poly_file_name_prefix)
  % execute the util triangle (http://www.cs.cmu.edu/~quake/triangle.html)
  % with the given command line arguements
  %
  % Input:
  % command_line_args: string of command line arguments
  %
  %
  % See corresponding paper: "Mixed finite elements for variational surface
  % modeling" by Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis Zorin, SGP
  % 2010
  %
  % Copyright 2010, Alec Jacobson, NYU
  %

  % hard code this for my computer for now, need to find out how to find
  % triangle...'which triangle' does not work
  path_to_triangle = '/opt/local/bin/triangle';

  [status, result] =system( ...
    [path_to_triangle ' ' command_line_args ' ' poly_file_name_prefix '.poly']);
  %status
  %result
  F = read_faces_from_ele_file([poly_file_name_prefix '.1.ele']);
  V = read_vertices_from_node_file([poly_file_name_prefix '.1.node']);
end
