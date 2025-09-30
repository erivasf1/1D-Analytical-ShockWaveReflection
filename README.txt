The problem set-up consists as follows:
  1) a stationary,compressible fluid w/ id of 1, originating at the left-end of the domain
  2) a stationary, compressible fluid w/ id of 2 that is to the right of fluid id 1 in which a moving piston separates the two fluids
  3) a solid wall at the right-end of the domain

The current assumptions are:
  1) fluid is a perfect gas
  2) piston moves at a constant velocity throughout the whole computation until it makes contact with the solid wall
  3) conservation of mass, momentum, and energy is satisfied across the shock (fluid 2), therefore Rankine-Hugoniot Jump conditions are employed to solve for intermediate states
  4) isentropic conditions are satisfied throughout the rarefaction fan (fluid 1)
