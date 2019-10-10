! Pre-multiply x weight by variable for speed
gx(-2) = gx(-2) * scaled_part_var
gx(-1) = gx(-1) * scaled_part_var
gx( 0) = gx( 0) * scaled_part_var
gx( 1) = gx( 1) * scaled_part_var
gx( 2) = gx( 2) * scaled_part_var

grid(cell_x1-2,cell_y1-2,cell_z1-2) = grid(cell_x1-2,cell_y1-2,cell_z1-2) &
    + gzy(-2,-2)*gx(-2)
grid(cell_x1-1,cell_y1-2,cell_z1-2) = grid(cell_x1-1,cell_y1-2,cell_z1-2) &
    + gzy(-2,-2)*gx(-1)
grid(cell_x1  ,cell_y1-2,cell_z1-2) = grid(cell_x1  ,cell_y1-2,cell_z1-2) &
    + gzy(-2,-2)*gx( 0)
grid(cell_x1+1,cell_y1-2,cell_z1-2) = grid(cell_x1+1,cell_y1-2,cell_z1-2) &
    + gzy(-2,-2)*gx( 1)
grid(cell_x1+2,cell_y1-2,cell_z1-2) = grid(cell_x1+2,cell_y1-2,cell_z1-2) &
    + gzy(-2,-2)*gx( 2)
grid(cell_x1-2,cell_y1-1,cell_z1-2) = grid(cell_x1-2,cell_y1-1,cell_z1-2) &
    + gzy(-2,-1)*gx(-2)
grid(cell_x1-1,cell_y1-1,cell_z1-2) = grid(cell_x1-1,cell_y1-1,cell_z1-2) &
    + gzy(-2,-1)*gx(-1)
grid(cell_x1  ,cell_y1-1,cell_z1-2) = grid(cell_x1  ,cell_y1-1,cell_z1-2) &
    + gzy(-2,-1)*gx( 0)
grid(cell_x1+1,cell_y1-1,cell_z1-2) = grid(cell_x1+1,cell_y1-1,cell_z1-2) &
    + gzy(-2,-1)*gx( 1)
grid(cell_x1+2,cell_y1-1,cell_z1-2) = grid(cell_x1+2,cell_y1-1,cell_z1-2) &
    + gzy(-2,-1)*gx( 2)
grid(cell_x1-2,cell_y1  ,cell_z1-2) = grid(cell_x1-2,cell_y1  ,cell_z1-2) &
    + gzy(-2, 0)*gx(-2)
grid(cell_x1-1,cell_y1  ,cell_z1-2) = grid(cell_x1-1,cell_y1  ,cell_z1-2) &
    + gzy(-2, 0)*gx(-1)
grid(cell_x1  ,cell_y1  ,cell_z1-2) = grid(cell_x1  ,cell_y1  ,cell_z1-2) &
    + gzy(-2, 0)*gx( 0)
grid(cell_x1+1,cell_y1  ,cell_z1-2) = grid(cell_x1+1,cell_y1  ,cell_z1-2) &
    + gzy(-2, 0)*gx( 1)
grid(cell_x1+2,cell_y1  ,cell_z1-2) = grid(cell_x1+2,cell_y1  ,cell_z1-2) &
    + gzy(-2, 0)*gx( 2)
grid(cell_x1-2,cell_y1+1,cell_z1-2) = grid(cell_x1-2,cell_y1+1,cell_z1-2) &
    + gzy(-2, 1)*gx(-2)
grid(cell_x1-1,cell_y1+1,cell_z1-2) = grid(cell_x1-1,cell_y1+1,cell_z1-2) &
    + gzy(-2, 1)*gx(-1)
grid(cell_x1  ,cell_y1+1,cell_z1-2) = grid(cell_x1  ,cell_y1+1,cell_z1-2) &
    + gzy(-2, 1)*gx( 0)
grid(cell_x1+1,cell_y1+1,cell_z1-2) = grid(cell_x1+1,cell_y1+1,cell_z1-2) &
    + gzy(-2, 1)*gx( 1)
grid(cell_x1+2,cell_y1+1,cell_z1-2) = grid(cell_x1+2,cell_y1+1,cell_z1-2) &
    + gzy(-2, 1)*gx( 2)
grid(cell_x1-2,cell_y1+2,cell_z1-2) = grid(cell_x1-2,cell_y1+2,cell_z1-2) &
    + gzy(-2, 2)*gx(-2)
grid(cell_x1-1,cell_y1+2,cell_z1-2) = grid(cell_x1-1,cell_y1+2,cell_z1-2) &
    + gzy(-2, 2)*gx(-1)
grid(cell_x1  ,cell_y1+2,cell_z1-2) = grid(cell_x1  ,cell_y1+2,cell_z1-2) &
    + gzy(-2, 2)*gx( 0)
grid(cell_x1+1,cell_y1+2,cell_z1-2) = grid(cell_x1+1,cell_y1+2,cell_z1-2) &
    + gzy(-2, 2)*gx( 1)
grid(cell_x1+2,cell_y1+2,cell_z1-2) = grid(cell_x1+2,cell_y1+2,cell_z1-2) &
    + gzy(-2, 2)*gx( 2)
grid(cell_x1-2,cell_y1-2,cell_z1-1) = grid(cell_x1-2,cell_y1-2,cell_z1-1) &
    + gzy(-1,-2)*gx(-2)
grid(cell_x1-1,cell_y1-2,cell_z1-1) = grid(cell_x1-1,cell_y1-2,cell_z1-1) &
    + gzy(-1,-2)*gx(-1)
grid(cell_x1  ,cell_y1-2,cell_z1-1) = grid(cell_x1  ,cell_y1-2,cell_z1-1) &
    + gzy(-1,-2)*gx( 0)
grid(cell_x1+1,cell_y1-2,cell_z1-1) = grid(cell_x1+1,cell_y1-2,cell_z1-1) &
    + gzy(-1,-2)*gx( 1)
grid(cell_x1+2,cell_y1-2,cell_z1-1) = grid(cell_x1+2,cell_y1-2,cell_z1-1) &
    + gzy(-1,-2)*gx( 2)
grid(cell_x1-2,cell_y1-1,cell_z1-1) = grid(cell_x1-2,cell_y1-1,cell_z1-1) &
    + gzy(-1,-1)*gx(-2)
grid(cell_x1-1,cell_y1-1,cell_z1-1) = grid(cell_x1-1,cell_y1-1,cell_z1-1) &
    + gzy(-1,-1)*gx(-1)
grid(cell_x1  ,cell_y1-1,cell_z1-1) = grid(cell_x1  ,cell_y1-1,cell_z1-1) &
    + gzy(-1,-1)*gx( 0)
grid(cell_x1+1,cell_y1-1,cell_z1-1) = grid(cell_x1+1,cell_y1-1,cell_z1-1) &
    + gzy(-1,-1)*gx( 1)
grid(cell_x1+2,cell_y1-1,cell_z1-1) = grid(cell_x1+2,cell_y1-1,cell_z1-1) &
    + gzy(-1,-1)*gx( 2)
grid(cell_x1-2,cell_y1  ,cell_z1-1) = grid(cell_x1-2,cell_y1  ,cell_z1-1) &
    + gzy(-1, 0)*gx(-2)
grid(cell_x1-1,cell_y1  ,cell_z1-1) = grid(cell_x1-1,cell_y1  ,cell_z1-1) &
    + gzy(-1, 0)*gx(-1)
grid(cell_x1  ,cell_y1  ,cell_z1-1) = grid(cell_x1  ,cell_y1  ,cell_z1-1) &
    + gzy(-1, 0)*gx( 0)
grid(cell_x1+1,cell_y1  ,cell_z1-1) = grid(cell_x1+1,cell_y1  ,cell_z1-1) &
    + gzy(-1, 0)*gx( 1)
grid(cell_x1+2,cell_y1  ,cell_z1-1) = grid(cell_x1+2,cell_y1  ,cell_z1-1) &
    + gzy(-1, 0)*gx( 2)
grid(cell_x1-2,cell_y1+1,cell_z1-1) = grid(cell_x1-2,cell_y1+1,cell_z1-1) &
    + gzy(-1, 1)*gx(-2)
grid(cell_x1-1,cell_y1+1,cell_z1-1) = grid(cell_x1-1,cell_y1+1,cell_z1-1) &
    + gzy(-1, 1)*gx(-1)
grid(cell_x1  ,cell_y1+1,cell_z1-1) = grid(cell_x1  ,cell_y1+1,cell_z1-1) &
    + gzy(-1, 1)*gx( 0)
grid(cell_x1+1,cell_y1+1,cell_z1-1) = grid(cell_x1+1,cell_y1+1,cell_z1-1) &
    + gzy(-1, 1)*gx( 1)
grid(cell_x1+2,cell_y1+1,cell_z1-1) = grid(cell_x1+2,cell_y1+1,cell_z1-1) &
    + gzy(-1, 1)*gx( 2)
grid(cell_x1-2,cell_y1+2,cell_z1-1) = grid(cell_x1-2,cell_y1+2,cell_z1-1) &
    + gzy(-1, 2)*gx(-2)
grid(cell_x1-1,cell_y1+2,cell_z1-1) = grid(cell_x1-1,cell_y1+2,cell_z1-1) &
    + gzy(-1, 2)*gx(-1)
grid(cell_x1  ,cell_y1+2,cell_z1-1) = grid(cell_x1  ,cell_y1+2,cell_z1-1) &
    + gzy(-1, 2)*gx( 0)
grid(cell_x1+1,cell_y1+2,cell_z1-1) = grid(cell_x1+1,cell_y1+2,cell_z1-1) &
    + gzy(-1, 2)*gx( 1)
grid(cell_x1+2,cell_y1+2,cell_z1-1) = grid(cell_x1+2,cell_y1+2,cell_z1-1) &
    + gzy(-1, 2)*gx( 2)
grid(cell_x1-2,cell_y1-2,cell_z1  ) = grid(cell_x1-2,cell_y1-2,cell_z1  ) &
    + gzy( 0,-2)*gx(-2)
grid(cell_x1-1,cell_y1-2,cell_z1  ) = grid(cell_x1-1,cell_y1-2,cell_z1  ) &
    + gzy( 0,-2)*gx(-1)
grid(cell_x1  ,cell_y1-2,cell_z1  ) = grid(cell_x1  ,cell_y1-2,cell_z1  ) &
    + gzy( 0,-2)*gx( 0)
grid(cell_x1+1,cell_y1-2,cell_z1  ) = grid(cell_x1+1,cell_y1-2,cell_z1  ) &
    + gzy( 0,-2)*gx( 1)
grid(cell_x1+2,cell_y1-2,cell_z1  ) = grid(cell_x1+2,cell_y1-2,cell_z1  ) &
    + gzy( 0,-2)*gx( 2)
grid(cell_x1-2,cell_y1-1,cell_z1  ) = grid(cell_x1-2,cell_y1-1,cell_z1  ) &
    + gzy( 0,-1)*gx(-2)
grid(cell_x1-1,cell_y1-1,cell_z1  ) = grid(cell_x1-1,cell_y1-1,cell_z1  ) &
    + gzy( 0,-1)*gx(-1)
grid(cell_x1  ,cell_y1-1,cell_z1  ) = grid(cell_x1  ,cell_y1-1,cell_z1  ) &
    + gzy( 0,-1)*gx( 0)
grid(cell_x1+1,cell_y1-1,cell_z1  ) = grid(cell_x1+1,cell_y1-1,cell_z1  ) &
    + gzy( 0,-1)*gx( 1)
grid(cell_x1+2,cell_y1-1,cell_z1  ) = grid(cell_x1+2,cell_y1-1,cell_z1  ) &
    + gzy( 0,-1)*gx( 2)
grid(cell_x1-2,cell_y1  ,cell_z1  ) = grid(cell_x1-2,cell_y1  ,cell_z1  ) &
    + gzy( 0, 0)*gx(-2)
grid(cell_x1-1,cell_y1  ,cell_z1  ) = grid(cell_x1-1,cell_y1  ,cell_z1  ) &
    + gzy( 0, 0)*gx(-1)
grid(cell_x1  ,cell_y1  ,cell_z1  ) = grid(cell_x1  ,cell_y1  ,cell_z1  ) &
    + gzy( 0, 0)*gx( 0)
grid(cell_x1+1,cell_y1  ,cell_z1  ) = grid(cell_x1+1,cell_y1  ,cell_z1  ) &
    + gzy( 0, 0)*gx( 1)
grid(cell_x1+2,cell_y1  ,cell_z1  ) = grid(cell_x1+2,cell_y1  ,cell_z1  ) &
    + gzy( 0, 0)*gx( 2)
grid(cell_x1-2,cell_y1+1,cell_z1  ) = grid(cell_x1-2,cell_y1+1,cell_z1  ) &
    + gzy( 0, 1)*gx(-2)
grid(cell_x1-1,cell_y1+1,cell_z1  ) = grid(cell_x1-1,cell_y1+1,cell_z1  ) &
    + gzy( 0, 1)*gx(-1)
grid(cell_x1  ,cell_y1+1,cell_z1  ) = grid(cell_x1  ,cell_y1+1,cell_z1  ) &
    + gzy( 0, 1)*gx( 0)
grid(cell_x1+1,cell_y1+1,cell_z1  ) = grid(cell_x1+1,cell_y1+1,cell_z1  ) &
    + gzy( 0, 1)*gx( 1)
grid(cell_x1+2,cell_y1+1,cell_z1  ) = grid(cell_x1+2,cell_y1+1,cell_z1  ) &
    + gzy( 0, 1)*gx( 2)
grid(cell_x1-2,cell_y1+2,cell_z1  ) = grid(cell_x1-2,cell_y1+2,cell_z1  ) &
    + gzy( 0, 2)*gx(-2)
grid(cell_x1-1,cell_y1+2,cell_z1  ) = grid(cell_x1-1,cell_y1+2,cell_z1  ) &
    + gzy( 0, 2)*gx(-1)
grid(cell_x1  ,cell_y1+2,cell_z1  ) = grid(cell_x1  ,cell_y1+2,cell_z1  ) &
    + gzy( 0, 2)*gx( 0)
grid(cell_x1+1,cell_y1+2,cell_z1  ) = grid(cell_x1+1,cell_y1+2,cell_z1  ) &
    + gzy( 0, 2)*gx( 1)
grid(cell_x1+2,cell_y1+2,cell_z1  ) = grid(cell_x1+2,cell_y1+2,cell_z1  ) &
    + gzy( 0, 2)*gx( 2)
grid(cell_x1-2,cell_y1-2,cell_z1+1) = grid(cell_x1-2,cell_y1-2,cell_z1+1) &
    + gzy( 1,-2)*gx(-2)
grid(cell_x1-1,cell_y1-2,cell_z1+1) = grid(cell_x1-1,cell_y1-2,cell_z1+1) &
    + gzy( 1,-2)*gx(-1)
grid(cell_x1  ,cell_y1-2,cell_z1+1) = grid(cell_x1  ,cell_y1-2,cell_z1+1) &
    + gzy( 1,-2)*gx( 0)
grid(cell_x1+1,cell_y1-2,cell_z1+1) = grid(cell_x1+1,cell_y1-2,cell_z1+1) &
    + gzy( 1,-2)*gx( 1)
grid(cell_x1+2,cell_y1-2,cell_z1+1) = grid(cell_x1+2,cell_y1-2,cell_z1+1) &
    + gzy( 1,-2)*gx( 2)
grid(cell_x1-2,cell_y1-1,cell_z1+1) = grid(cell_x1-2,cell_y1-1,cell_z1+1) &
    + gzy( 1,-1)*gx(-2)
grid(cell_x1-1,cell_y1-1,cell_z1+1) = grid(cell_x1-1,cell_y1-1,cell_z1+1) &
    + gzy( 1,-1)*gx(-1)
grid(cell_x1  ,cell_y1-1,cell_z1+1) = grid(cell_x1  ,cell_y1-1,cell_z1+1) &
    + gzy( 1,-1)*gx( 0)
grid(cell_x1+1,cell_y1-1,cell_z1+1) = grid(cell_x1+1,cell_y1-1,cell_z1+1) &
    + gzy( 1,-1)*gx( 1)
grid(cell_x1+2,cell_y1-1,cell_z1+1) = grid(cell_x1+2,cell_y1-1,cell_z1+1) &
    + gzy( 1,-1)*gx( 2)
grid(cell_x1-2,cell_y1  ,cell_z1+1) = grid(cell_x1-2,cell_y1  ,cell_z1+1) &
    + gzy( 1, 0)*gx(-2)
grid(cell_x1-1,cell_y1  ,cell_z1+1) = grid(cell_x1-1,cell_y1  ,cell_z1+1) &
    + gzy( 1, 0)*gx(-1)
grid(cell_x1  ,cell_y1  ,cell_z1+1) = grid(cell_x1  ,cell_y1  ,cell_z1+1) &
    + gzy( 1, 0)*gx( 0)
grid(cell_x1+1,cell_y1  ,cell_z1+1) = grid(cell_x1+1,cell_y1  ,cell_z1+1) &
    + gzy( 1, 0)*gx( 1)
grid(cell_x1+2,cell_y1  ,cell_z1+1) = grid(cell_x1+2,cell_y1  ,cell_z1+1) &
    + gzy( 1, 0)*gx( 2)
grid(cell_x1-2,cell_y1+1,cell_z1+1) = grid(cell_x1-2,cell_y1+1,cell_z1+1) &
    + gzy( 1, 1)*gx(-2)
grid(cell_x1-1,cell_y1+1,cell_z1+1) = grid(cell_x1-1,cell_y1+1,cell_z1+1) &
    + gzy( 1, 1)*gx(-1)
grid(cell_x1  ,cell_y1+1,cell_z1+1) = grid(cell_x1  ,cell_y1+1,cell_z1+1) &
    + gzy( 1, 1)*gx( 0)
grid(cell_x1+1,cell_y1+1,cell_z1+1) = grid(cell_x1+1,cell_y1+1,cell_z1+1) &
    + gzy( 1, 1)*gx( 1)
grid(cell_x1+2,cell_y1+1,cell_z1+1) = grid(cell_x1+2,cell_y1+1,cell_z1+1) &
    + gzy( 1, 1)*gx( 2)
grid(cell_x1-2,cell_y1+2,cell_z1+1) = grid(cell_x1-2,cell_y1+2,cell_z1+1) &
    + gzy( 1, 2)*gx(-2)
grid(cell_x1-1,cell_y1+2,cell_z1+1) = grid(cell_x1-1,cell_y1+2,cell_z1+1) &
    + gzy( 1, 2)*gx(-1)
grid(cell_x1  ,cell_y1+2,cell_z1+1) = grid(cell_x1  ,cell_y1+2,cell_z1+1) &
    + gzy( 1, 2)*gx( 0)
grid(cell_x1+1,cell_y1+2,cell_z1+1) = grid(cell_x1+1,cell_y1+2,cell_z1+1) &
    + gzy( 1, 2)*gx( 1)
grid(cell_x1+2,cell_y1+2,cell_z1+1) = grid(cell_x1+2,cell_y1+2,cell_z1+1) &
    + gzy( 1, 2)*gx( 2)
grid(cell_x1-2,cell_y1-2,cell_z1+2) = grid(cell_x1-2,cell_y1-2,cell_z1+2) &
    + gzy( 2,-2)*gx(-2)
grid(cell_x1-1,cell_y1-2,cell_z1+2) = grid(cell_x1-1,cell_y1-2,cell_z1+2) &
    + gzy( 2,-2)*gx(-1)
grid(cell_x1  ,cell_y1-2,cell_z1+2) = grid(cell_x1  ,cell_y1-2,cell_z1+2) &
    + gzy( 2,-2)*gx( 0)
grid(cell_x1+1,cell_y1-2,cell_z1+2) = grid(cell_x1+1,cell_y1-2,cell_z1+2) &
    + gzy( 2,-2)*gx( 1)
grid(cell_x1+2,cell_y1-2,cell_z1+2) = grid(cell_x1+2,cell_y1-2,cell_z1+2) &
    + gzy( 2,-2)*gx( 2)
grid(cell_x1-2,cell_y1-1,cell_z1+2) = grid(cell_x1-2,cell_y1-1,cell_z1+2) &
    + gzy( 2,-1)*gx(-2)
grid(cell_x1-1,cell_y1-1,cell_z1+2) = grid(cell_x1-1,cell_y1-1,cell_z1+2) &
    + gzy( 2,-1)*gx(-1)
grid(cell_x1  ,cell_y1-1,cell_z1+2) = grid(cell_x1  ,cell_y1-1,cell_z1+2) &
    + gzy( 2,-1)*gx( 0)
grid(cell_x1+1,cell_y1-1,cell_z1+2) = grid(cell_x1+1,cell_y1-1,cell_z1+2) &
    + gzy( 2,-1)*gx( 1)
grid(cell_x1+2,cell_y1-1,cell_z1+2) = grid(cell_x1+2,cell_y1-1,cell_z1+2) &
    + gzy( 2,-1)*gx( 2)
grid(cell_x1-2,cell_y1  ,cell_z1+2) = grid(cell_x1-2,cell_y1  ,cell_z1+2) &
    + gzy( 2, 0)*gx(-2)
grid(cell_x1-1,cell_y1  ,cell_z1+2) = grid(cell_x1-1,cell_y1  ,cell_z1+2) &
    + gzy( 2, 0)*gx(-1)
grid(cell_x1  ,cell_y1  ,cell_z1+2) = grid(cell_x1  ,cell_y1  ,cell_z1+2) &
    + gzy( 2, 0)*gx( 0)
grid(cell_x1+1,cell_y1  ,cell_z1+2) = grid(cell_x1+1,cell_y1  ,cell_z1+2) &
    + gzy( 2, 0)*gx( 1)
grid(cell_x1+2,cell_y1  ,cell_z1+2) = grid(cell_x1+2,cell_y1  ,cell_z1+2) &
    + gzy( 2, 0)*gx( 2)
grid(cell_x1-2,cell_y1+1,cell_z1+2) = grid(cell_x1-2,cell_y1+1,cell_z1+2) &
    + gzy( 2, 1)*gx(-2)
grid(cell_x1-1,cell_y1+1,cell_z1+2) = grid(cell_x1-1,cell_y1+1,cell_z1+2) &
    + gzy( 2, 1)*gx(-1)
grid(cell_x1  ,cell_y1+1,cell_z1+2) = grid(cell_x1  ,cell_y1+1,cell_z1+2) &
    + gzy( 2, 1)*gx( 0)
grid(cell_x1+1,cell_y1+1,cell_z1+2) = grid(cell_x1+1,cell_y1+1,cell_z1+2) &
    + gzy( 2, 1)*gx( 1)
grid(cell_x1+2,cell_y1+1,cell_z1+2) = grid(cell_x1+2,cell_y1+1,cell_z1+2) &
    + gzy( 2, 1)*gx( 2)
grid(cell_x1-2,cell_y1+2,cell_z1+2) = grid(cell_x1-2,cell_y1+2,cell_z1+2) &
    + gzy( 2, 2)*gx(-2)
grid(cell_x1-1,cell_y1+2,cell_z1+2) = grid(cell_x1-1,cell_y1+2,cell_z1+2) &
    + gzy( 2, 2)*gx(-1)
grid(cell_x1  ,cell_y1+2,cell_z1+2) = grid(cell_x1  ,cell_y1+2,cell_z1+2) &
    + gzy( 2, 2)*gx( 0)
grid(cell_x1+1,cell_y1+2,cell_z1+2) = grid(cell_x1+1,cell_y1+2,cell_z1+2) &
    + gzy( 2, 2)*gx( 1)
grid(cell_x1+2,cell_y1+2,cell_z1+2) = grid(cell_x1+2,cell_y1+2,cell_z1+2) &
    + gzy( 2, 2)*gx( 2)
