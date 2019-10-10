! Pre-multiply x weight by variable for speed
gx(-1) = gx(-1) * scaled_part_var
gx( 0) = gx( 0) * scaled_part_var
gx( 1) = gx( 1) * scaled_part_var

grid(cell_x1-1,cell_y1-1,cell_z1-1) = grid(cell_x1-1,cell_y1-1,cell_z1-1) &
    + gzy(-1,-1) * gx(-1)
grid(cell_x1  ,cell_y1-1,cell_z1-1) = grid(cell_x1  ,cell_y1-1,cell_z1-1) &
    + gzy(-1,-1) * gx( 0)
grid(cell_x1+1,cell_y1-1,cell_z1-1) = grid(cell_x1+1,cell_y1-1,cell_z1-1) &
    + gzy(-1,-1) * gx( 1)
grid(cell_x1-1,cell_y1  ,cell_z1-1) = grid(cell_x1-1,cell_y1  ,cell_z1-1) &
    + gzy(-1, 0) * gx(-1)
grid(cell_x1  ,cell_y1  ,cell_z1-1) = grid(cell_x1  ,cell_y1  ,cell_z1-1) &
    + gzy(-1, 0) * gx( 0)
grid(cell_x1+1,cell_y1  ,cell_z1-1) = grid(cell_x1+1,cell_y1  ,cell_z1-1) &
    + gzy(-1, 0) * gx( 1)
grid(cell_x1-1,cell_y1+1,cell_z1-1) = grid(cell_x1-1,cell_y1+1,cell_z1-1) &
    + gzy(-1, 1) * gx(-1)
grid(cell_x1  ,cell_y1+1,cell_z1-1) = grid(cell_x1  ,cell_y1+1,cell_z1-1) &
    + gzy(-1, 1) * gx( 0)
grid(cell_x1+1,cell_y1+1,cell_z1-1) = grid(cell_x1+1,cell_y1+1,cell_z1-1) &
    + gzy(-1, 1) * gx( 1)
grid(cell_x1-1,cell_y1-1,cell_z1  ) = grid(cell_x1-1,cell_y1-1,cell_z1  ) &
    + gzy( 0,-1) * gx(-1)
grid(cell_x1  ,cell_y1-1,cell_z1  ) = grid(cell_x1  ,cell_y1-1,cell_z1  ) &
    + gzy( 0,-1) * gx( 0)
grid(cell_x1+1,cell_y1-1,cell_z1  ) = grid(cell_x1+1,cell_y1-1,cell_z1  ) &
    + gzy( 0,-1) * gx( 1)
grid(cell_x1-1,cell_y1  ,cell_z1  ) = grid(cell_x1-1,cell_y1  ,cell_z1  ) &
    + gzy( 0, 0) * gx(-1)
grid(cell_x1  ,cell_y1  ,cell_z1  ) = grid(cell_x1  ,cell_y1  ,cell_z1  ) &
    + gzy( 0, 0) * gx( 0)
grid(cell_x1+1,cell_y1  ,cell_z1  ) = grid(cell_x1+1,cell_y1  ,cell_z1  ) &
    + gzy( 0, 0) * gx( 1)
grid(cell_x1-1,cell_y1+1,cell_z1  ) = grid(cell_x1-1,cell_y1+1,cell_z1  ) &
    + gzy( 0, 1) * gx(-1)
grid(cell_x1  ,cell_y1+1,cell_z1  ) = grid(cell_x1  ,cell_y1+1,cell_z1  ) &
    + gzy( 0, 1) * gx( 0)
grid(cell_x1+1,cell_y1+1,cell_z1  ) = grid(cell_x1+1,cell_y1+1,cell_z1  ) &
    + gzy( 0, 1) * gx( 1)
grid(cell_x1-1,cell_y1-1,cell_z1+1) = grid(cell_x1-1,cell_y1-1,cell_z1+1) &
    + gzy( 1,-1) * gx(-1)
grid(cell_x1  ,cell_y1-1,cell_z1+1) = grid(cell_x1  ,cell_y1-1,cell_z1+1) &
    + gzy( 1,-1) * gx( 0)
grid(cell_x1+1,cell_y1-1,cell_z1+1) = grid(cell_x1+1,cell_y1-1,cell_z1+1) &
    + gzy( 1,-1) * gx( 1)
grid(cell_x1-1,cell_y1  ,cell_z1+1) = grid(cell_x1-1,cell_y1  ,cell_z1+1) &
    + gzy( 1, 0) * gx(-1)
grid(cell_x1  ,cell_y1  ,cell_z1+1) = grid(cell_x1  ,cell_y1  ,cell_z1+1) &
    + gzy( 1, 0) * gx( 0)
grid(cell_x1+1,cell_y1  ,cell_z1+1) = grid(cell_x1+1,cell_y1  ,cell_z1+1) &
    + gzy( 1, 0) * gx( 1)
grid(cell_x1-1,cell_y1+1,cell_z1+1) = grid(cell_x1-1,cell_y1+1,cell_z1+1) &
    + gzy( 1, 1) * gx(-1)
grid(cell_x1  ,cell_y1+1,cell_z1+1) = grid(cell_x1  ,cell_y1+1,cell_z1+1) &
    + gzy( 1, 1) * gx( 0)
grid(cell_x1+1,cell_y1+1,cell_z1+1) = grid(cell_x1+1,cell_y1+1,cell_z1+1) &
    + gzy( 1, 1) * gx( 1)
