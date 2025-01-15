
function EnergyBalanceEquation(x) result(fx)
    !Define the energy equation to be solved by Brent's method
    implicit none
    real :: x, fx !x is the root to be determined, fx is the value to be optimized.
    
    !fx = c* x**2 + d*x +e ! explicit expression of the scalar equation
    !x=TZERO(I)
    !fx=RESID(I) 


    WZERO(I) = 0.622 * calcEsat(x) / PADRY(I)
    Q0SAT(I) = WZERO(I) / (1.0 + WZERO(I))
    if (IWATER(I) > 0) then
      EVBETA(I) = 1.0
      QZERO(I) = Q0SAT(I)
    else
      if (ipeatland(i) == 0) then
        EVBETA(I) = CEVAP(I)
      else
        ievap(i) = ievapmoss(i)
        evbeta(i) = cevapmoss(i)
      end if
      QZERO(I) = EVBETA(I) * Q0SAT(I) + (1.0 - EVBETA(I)) * QAC(I)
      if (QZERO(I) > QAC(I) .and. IEVAP(I) == 0) then
        EVBETA(I) = 0.0
        QZERO(I) = QAC(I)
      end if
    end if
    !
    TPOTG(I) = x - 8.0 * ZOM(I) * GRAV / SPHAIR
    TVIRTG(I) = TPOTG(I) * (1.0 + 0.61 * QZERO(I))
    if (TVIRTG(I) > TVRTAC(I) + 1.) then
      RAGINV(I) = RAGCO * (TVIRTG(I) - TVRTAC(I)) ** 0.333333
      DRAGIN(I) = 0.333 * RAGCO * (TVIRTG(I) - TVRTAC(I)) ** ( - .667)
    else if (TVIRTG(I) > (TVRTAC(I) + 0.001)) then
      RAGINV(I) = RAGCO * (TVIRTG(I) - TVRTAC(I))
      DRAGIN(I) = RAGCO
    else
      RAGINV(I) = 0.0
      DRAGIN(I) = 0.0
    end if
    !
    QLWOG(I) = SBC * x * x * x * x
    QSENSG(I) = RHOAIR(I) * SPHAIR * RAGINV(I) * &
                (TPOTG(I) - TAC(I))
    EVAPG (I) = RHOAIR(I) * (QZERO(I) - QAC(I)) * RAGINV(I)
    if (EVAPG(I) > EVPMAX(I)) EVAPG(I) = EVPMAX(I)
    QEVAPG(I) = CPHCHG(I) * EVAPG(I)
    GZERO(I) = GCOEFF(I) * x + GCONST(I)
    fx = QSWNG(I) + FSVF(I) * QLWIN(I) + (1.0 - FSVF(I)) * &
                QLWOC(I) - QLWOG(I) - QSENSG(I) - QEVAPG(I) - GZERO(I)

  end function EnergyBalanceEquation