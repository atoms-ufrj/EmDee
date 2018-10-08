block
  real(rb) :: r2fac, u, u2, u3, G, WG

  select case (model%modifier)
#   if defined(compute)
      case (SHIFTED)
        Eij = Eij + model%eshift
#   endif
    case (SHIFTED_FORCE)
      rFc = model%fshift/invR
      Wij = Wij - rFc
#     if defined(compute)
        Eij = Eij + model%eshift + rFc
#     endif
    case (SMOOTHED, SHIFTED_SMOOTHED)
      r2fac = model%factor/invR
      if (r2fac > model%Rm2fac) then
#       ifndef compute
          select type ( model )
            include energy_file
          end select
#       endif
        u = r2fac - model%Rm2fac
        u2 = u*u
        u3 = u*u2
        G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
        WG = -30.0_rb*u2*(2.0_rb*u - u2 - 1.0_rb)*r2fac
        Eij = Eij + model%eshift
        Wij = Wij*G + Eij*WG
        Eij = Eij*G
      end if
    case (SQUARE_SMOOTHED, SHIFTED_SQUARE_SMOOTHED)
      r2fac = model%factor/invR2
      if (r2fac > model%Rm2fac) then
#       ifndef compute
          select type ( model )
            include energy_file
          end select
#       endif
        u = r2fac - model%Rm2fac
        u2 = u*u
        u3 = u*u2
        G = 1.0_rb + u3*(15.0_rb*u - 6.0_rb*u2 - 10.0_rb)
        WG = -60.0_rb*u2*(2.0_rb*u - u2 - 1.0_rb)*r2fac
        Eij = Eij + model%eshift
        Wij = Wij*G + Eij*WG
        Eij = Eij*G
      end if
  end select
end block
