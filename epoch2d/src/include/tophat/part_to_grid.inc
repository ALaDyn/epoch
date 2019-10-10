! Pre-multiply x weight by variable for speed
gx( 0) = gx( 0) * scaled_part_var
gx( 1) = gx( 1) * scaled_part_var

grid(cell_x1  ,cell_y1  ) = grid(cell_x1  ,cell_y1  ) + gy(0) * gx(0)
grid(cell_x1+1,cell_y1  ) = grid(cell_x1+1,cell_y1  ) + gy(0) * gx(1)
grid(cell_x1  ,cell_y1+1) = grid(cell_x1  ,cell_y1+1) + gy(1) * gx(0)
grid(cell_x1+1,cell_y1+1) = grid(cell_x1+1,cell_y1+1) + gy(1) * gx(1)
