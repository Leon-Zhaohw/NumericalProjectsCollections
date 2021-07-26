function [ nodes ] = cellsToNodes( cells )
%CELLSTONODES Convert cells to nodes. Both are logical matrices. A single 
%       true cell incident to the node suffices to set the node to true. 
%       This prevents fluids and solids from being only a single node wide.
%   return nodes boolean matrix(yCells + 1, xCells + 1)
%   cells: boolean matrix(yCells, xCells)

    yNodes = size(cells, 1) + 1;
    xNodes = size(cells, 2) + 1;
    
    % for easier combination of matrizes
    paddedCells = padarray(cells, [1 1], 'replicate');
    cellsOneUp = paddedCells(1 : yNodes, 2 : xNodes + 1);
    cellsOneCenter = paddedCells(1 : yNodes, 1 : xNodes);
    cellsOneRight = paddedCells(2 : yNodes + 1, 1 : xNodes);
    cellsOneRightUp = paddedCells(2 : yNodes + 1, 2 : xNodes + 1);
    
    nodes = cellsOneUp | cellsOneCenter | cellsOneRightUp | cellsOneRight;

end
